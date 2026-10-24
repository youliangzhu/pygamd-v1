'''
PYGAMD eager-graph execution policy resolver.

This module is intentionally CUDA-driver-free.  It centralizes the new primary
runtime policy names while preserving the legacy low-level environment knobs
used by the external validation scripts.
'''

import os

POLICY_ENV = "PYGAMD_EXECUTION_POLICY"
DEFAULT_POLICY = "eager_graph"
VALID_POLICIES = ("eager_graph", "graph_strict", "eager_debug")

# Execution route tiers: graph is the primary execution route; pure eager is retained
# only as an explicit debugging / reference baseline path.  These labels are
# intentionally surfaced in diagnostics and benchmark JSON so reports do not
# accidentally present pure eager as the official CuPy route.
POLICY_ROUTE_TIER = {
	"eager_graph": "primary_eager_graph",
	"graph_strict": "coverage_ci_graph",
	"eager_debug": "legacy_eager_debug",
}
POLICY_ROUTE_LABEL = {
	"eager_graph": "eager-graph primary route",
	"graph_strict": "graph-strict coverage/CI route",
	"eager_debug": "legacy pure-eager debug/reference route",
}

_FALSEY = ("", "0", "false", "off", "no", "n")
_POLICY_ALIASES = {
	"": DEFAULT_POLICY,
	"default": DEFAULT_POLICY,
	"auto": DEFAULT_POLICY,
	"primary": DEFAULT_POLICY,
	"graph": DEFAULT_POLICY,
	"eager_graph": "eager_graph",
	"eagergraph": "eager_graph",
	"fusion_graph": "eager_graph",
	"cuda_graph": "eager_graph",
	"graph_primary": "eager_graph",
	"strict": "graph_strict",
	"graph_strict": "graph_strict",
	"strict_graph": "graph_strict",
	"coverage": "graph_strict",
	"ci": "graph_strict",
	"debug": "eager_debug",
	"eager": "eager_debug",
	"pure_eager": "eager_debug",
	"legacy_eager": "eager_debug",
	"eager_debug": "eager_debug",
	"off": "eager_debug",
}


def _clean(value):
	return str(value if value is not None else "").strip().lower().replace("-", "_")


def _env_is_explicit(name):
	return name in os.environ and str(os.environ.get(name, "")).strip() != ""


def _legacy_cuda_graph_requests_off():
	if "PYGAMD_CUDA_GRAPH" not in os.environ:
		return False
	v = _clean(os.environ.get("PYGAMD_CUDA_GRAPH"))
	return v in _FALSEY


def _parse_float_env(name, default):
	raw = os.environ.get(name, "")
	if str(raw).strip() == "":
		return default
	try:
		value = float(raw)
	except Exception:
		return default
	if value > 1.0:
		value = value / 100.0
	if value < 0.0:
		return default
	if value > 1.0:
		return 1.0
	return value


def _parse_int_env(name, default=None):
	raw = os.environ.get(name, "")
	if str(raw).strip() == "":
		return default
	try:
		return int(raw)
	except Exception:
		return default


class ExecutionPolicy:
	'''Resolved pygamd runtime policy.

	name:
	  eager_graph  -> graph primary route, fallback/eager islands allowed.
	  graph_strict -> graph primary route plus runtime coverage guard.
	  eager_debug  -> legacy pure-eager debug route; graph/fusion disabled.
	'''
	__slots__ = ("name", "explicit", "graph_enabled", "graph_strict", "allow_eager_fallback", "route_tier", "route_label", "is_primary_route", "is_legacy_debug_route")

	def __init__(self, name, explicit=False):
		self.name = name
		self.explicit = bool(explicit)
		self.graph_enabled = name in ("eager_graph", "graph_strict")
		self.graph_strict = (name == "graph_strict")
		self.allow_eager_fallback = (name != "graph_strict")
		self.route_tier = POLICY_ROUTE_TIER.get(name, "unknown")
		self.route_label = POLICY_ROUTE_LABEL.get(name, name)
		self.is_primary_route = (name == "eager_graph")
		self.is_legacy_debug_route = (name == "eager_debug")

	def as_dict(self):
		return {
			"pygamd_execution_policy": self.name,
			"pygamd_execution_policy_explicit": self.explicit,
			"pygamd_execution_graph_enabled": self.graph_enabled,
			"pygamd_execution_graph_strict": self.graph_strict,
			"pygamd_execution_allow_eager_fallback": self.allow_eager_fallback,
			"pygamd_execution_route_tier": self.route_tier,
			"pygamd_execution_route_label": self.route_label,
			"pygamd_execution_is_primary_route": self.is_primary_route,
			"pygamd_execution_is_legacy_debug_route": self.is_legacy_debug_route,
		}

	def __repr__(self):
		return "ExecutionPolicy(name=%r, explicit=%r)" % (self.name, self.explicit)


def normalize_policy(value=None):
	'''Return one of VALID_POLICIES.

	When PYGAMD_EXECUTION_POLICY is unset, the new default is eager_graph.  The only
	legacy exception is an explicit PYGAMD_CUDA_GRAPH=0/off, which is interpreted as
	a request for the hidden eager_debug route so old debugging commands still work.
	'''
	explicit = value is not None or _env_is_explicit(POLICY_ENV)
	raw = os.environ.get(POLICY_ENV, "") if value is None else value
	key = _clean(raw)
	if not explicit and key == "":
		if _legacy_cuda_graph_requests_off():
			return "eager_debug"
		return DEFAULT_POLICY
	if key in _POLICY_ALIASES:
		return _POLICY_ALIASES[key]
	raise RuntimeError(
		"Unsupported %s=%r. Valid policies: eager_graph, graph_strict, eager_debug."
		% (POLICY_ENV, raw))


def resolve_execution_policy(value=None):
	return ExecutionPolicy(normalize_policy(value), explicit=(value is not None or _env_is_explicit(POLICY_ENV)))


def apply_execution_policy_defaults(policy=None):
	'''Apply policy-owned defaults for legacy low-level knobs.

	The graph controller uses this once at dynamics construction time.  The updates
	are deliberate: RawNVT fusion and DPD graph dynamic seed are still read by their
	own modules from legacy env vars, so the high-level policy must seed those knobs
	before the first kernel path is initialized.
	'''
	policy = policy or resolve_execution_policy()
	if policy.name == "eager_debug":
		os.environ["PYGAMD_CUDA_GRAPH"] = "0"
		# Default to the legacy pure-eager debug route, but do not clobber an
		# explicit low-level fusion/dynamic-seed knob used by old comparison scripts.
		if "PYGAMD_RAWNVT_FUSION" not in os.environ or str(os.environ.get("PYGAMD_RAWNVT_FUSION", "")).strip() == "":
			os.environ["PYGAMD_RAWNVT_FUSION"] = "0"
		if "PYGAMD_DPD_GRAPH_DYNAMIC_SEED" not in os.environ or str(os.environ.get("PYGAMD_DPD_GRAPH_DYNAMIC_SEED", "")).strip() == "":
			os.environ["PYGAMD_DPD_GRAPH_DYNAMIC_SEED"] = "0"
		return policy

	# eager_graph / graph_strict are the new primary graph-capable routes.
	# Keep an explicit legacy "tail" request as a mode override; otherwise ensure
	# graph is on even when old scripts forgot to export PYGAMD_CUDA_GRAPH=1.
	if _clean(os.environ.get("PYGAMD_CUDA_GRAPH", "")) != "tail":
		os.environ["PYGAMD_CUDA_GRAPH"] = "1"
	if "PYGAMD_RAWNVT_FUSION" not in os.environ or str(os.environ.get("PYGAMD_RAWNVT_FUSION", "")).strip() == "":
		os.environ["PYGAMD_RAWNVT_FUSION"] = "1"
	if "PYGAMD_DPD_GRAPH_DYNAMIC_SEED" not in os.environ or str(os.environ.get("PYGAMD_DPD_GRAPH_DYNAMIC_SEED", "")).strip() == "":
		os.environ["PYGAMD_DPD_GRAPH_DYNAMIC_SEED"] = "1"
	if policy.graph_strict and "PYGAMD_CUDA_GRAPH_DIAG" not in os.environ:
		os.environ["PYGAMD_CUDA_GRAPH_DIAG"] = "1"
	return policy


def policy_route_tier(policy=None):
	policy = policy or resolve_execution_policy()
	return getattr(policy, "route_tier", POLICY_ROUTE_TIER.get(getattr(policy, "name", policy), "unknown"))


def policy_route_label(policy=None):
	policy = policy or resolve_execution_policy()
	return getattr(policy, "route_label", POLICY_ROUTE_LABEL.get(getattr(policy, "name", policy), str(policy)))


def is_legacy_eager_debug_policy(policy=None):
	policy = policy or resolve_execution_policy()
	return bool(getattr(policy, "is_legacy_debug_route", False))


def read_cuda_graph_mode(policy=None):
	policy = policy or resolve_execution_policy()
	if not policy.graph_enabled:
		return False, "off"
	v = _clean(os.environ.get("PYGAMD_CUDA_GRAPH", "1"))
	if v == "tail":
		return True, "tail"
	if v in _FALSEY:
		# Explicit high-level graph policies win over stale legacy off values.
		if policy.explicit:
			return True, _read_cuda_graph_mode_value()
		return False, "off"
	return True, _read_cuda_graph_mode_value()


def _read_cuda_graph_mode_value():
	m = _clean(os.environ.get("PYGAMD_CUDA_GRAPH_MODE", "full"))
	return "tail" if m == "tail" else "full"


def graph_strict_min_coverage(default=0.90):
	return _parse_float_env("PYGAMD_GRAPH_STRICT_MIN_COVERAGE", default)


def graph_strict_max_eager_steps(default=None):
	return _parse_int_env("PYGAMD_GRAPH_STRICT_MAX_EAGER_STEPS", default)
