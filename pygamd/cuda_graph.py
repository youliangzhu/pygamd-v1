'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License.
	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE.
	See the General Public License v3 for more details.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''


# -----------------------------------------------------------------------------
# CUDA Graph stable-segment engine.
#
# Execution layers per step (application orchestrates):
#   1. FULL graph  : firststep + force segment + temp_partial + xi_update +
#                    secondstep -> one graph launch.  Eligible only on predicted
#                    nlist-skip steps.
#   2. TAIL graph  : eager firststep + eager nlist + graph(force segment + NVT).
#   3. EAGER       : original path.
#
# CUDA Graph broadens the graph force segment while preserving fallback safety:
#   * single RawPair still uses the old overwrite force graph (no zero-force fill);
#   * multiple graph-compatible additive forces use captured zero_force_vp +
#     sequential graph_force_add() launches;
#   * conservative DPD (sigma == 0 for all pairs) can be graph-captured;
#   * CUDA Graph stochastic DPD can opt in to stochastic/dissipative DPD graph replay by
#     updating device-side scalar inputs before each graph.launch(), preserving
#     Saru(base_seed+timestep, min/max tag) and variant temperature schedules
#     without freezing by-value scalars in the captured graph.
#
# Env: PYGAMD_EXECUTION_POLICY = eager_graph (default) | graph_strict | eager_debug
#      PYGAMD_CUDA_GRAPH remains as a low-level compatibility knob/mode override:
#        0/off -> eager_debug when no high-level policy is set; 1/auto -> full; tail -> tail
#      PYGAMD_CUDA_GRAPH_MODE = full (default when graph on) | tail
#      PYGAMD_GRAPH_STRICT_MIN_COVERAGE = 0.90 by default for graph_strict
#      PYGAMD_CUDA_GRAPH_PROFILE = 1   (print capture/replay/fallback summary)
#      PYGAMD_CUDA_GRAPH_DIAG = 1      (structured engine metrics diagnostics)
# Anything unsupported safely falls back to eager unless graph_strict coverage fails.
# -----------------------------------------------------------------------------

import os

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp
from pygamd import engine_timing as _t
from pygamd.execution_policy import (
	apply_execution_policy_defaults,
	graph_strict_max_eager_steps,
	graph_strict_min_coverage,
	read_cuda_graph_mode,
	resolve_execution_policy,
)
from pygamd.graph_protocol import (
	FORCE_GRAPH_PROTOCOL_VERSION,
	describe_force_graph_backend,
	force_graph_backend_from_wrapper,
	force_graph_backend_name,
	force_graph_can_add,
	force_graph_can_write,
	force_graph_prepare_replay,
	force_graph_ready,
	force_graph_signature,
	force_graph_unready_reason,
)


def _read_mode(policy=None):
	return read_cuda_graph_mode(policy)


def _gprofile_enabled():
	v = os.environ.get("PYGAMD_CUDA_GRAPH_PROFILE", "")
	return v not in ("", "0", "false", "False", "no", "off")


def _diag_enabled():
	v = os.environ.get("PYGAMD_CUDA_GRAPH_DIAG", "")
	return v not in ("", "0", "false", "False", "no", "off")


def _sp_inc(name, n=1):
	if _sp.ENABLED:
		_sp.inc(name, n)


class StableCudaGraph:
	def __init__(self):
		self.policy = apply_execution_policy_defaults(resolve_execution_policy())
		self.execution_policy = self.policy.name
		self.strict = bool(self.policy.graph_strict)
		self.strict_min_coverage = graph_strict_min_coverage()
		self.strict_max_eager_steps = graph_strict_max_eager_steps()
		self.enabled, self.mode = _read_mode(self.policy)   # mode: off | tail | full
		self.allow_full = (self.mode == "full")
		self.gprofile = _gprofile_enabled()
		self.gdiag = _diag_enabled() or self.strict
		self.cp = None
		self.available = None
		self.stream = None
		self._disabled = False
		self._tail_plan = None
		self.application_loop_kind = "not_started"
		self.application_loop_label = "not started"
		self.application_loop_count = 0
		self.application_loop_eager_graph_count = 0
		self.application_loop_eager_debug_count = 0
		# Legacy aliases retained for existing external result consumers.
		self.application_loop = self.application_loop_kind
		self.application_loop_graph_aware = False
		self.force_graph_protocol_version = FORCE_GRAPH_PROTOCOL_VERSION
		self.force_graph_protocol_summary = []

		# full graph
		self.full_graph = None
		self.full_sig = None
		self.full_target = None
		# tail graph
		self.tail_graph = None
		self.tail_sig = None
		self.tail_target = None

		# counters
		self.full_capture = 0
		self.full_replay = 0
		self.tail_capture = 0
		self.tail_replay = 0
		self.full_step_count = 0
		self.tail_step_count = 0
		self.eager_step_count = 0
		self.fallback_count = 0
		self.reasons = {}

	# ---- counters ----
	def _fb(self, name):
		name = str(name or "unsupported")
		self.reasons[name] = self.reasons.get(name, 0) + 1
		self.fallback_count += 1
		prof.count("cuda_graph_fallback_count")
		prof.count("cuda_graph_fallback_" + name)
		_sp_inc("cuda_graph_fallback_count")
		_sp_inc("cuda_graph_fallback_reason_" + name)

	def record_eager_step(self):
		'''Record a whole-step force/integrator-tail eager fallback.

		This deliberately does NOT count controlled eager islands such as nlist
		checks/rebuilds, dumps, sort, or the eager firststep used by tail graphs.
		graph_strict uses this count to detect cases where the graph route has been
		replaced by the legacy pure-eager force/tail path.
		'''
		if not self.enabled:
			return
		self.eager_step_count += 1
		prof.count("cuda_graph_eager_step")
		_sp_inc("cuda_graph_eager_step_count")

	def _prepare_graph_inputs(self, plan, timestep):
		"""Update graph-dynamic device inputs immediately before replay.

		Captured graphs hold stable device pointers for dynamic scalar inputs.
		Current users include BD/DPD random seeds and NVT/BD/DPD variant
		temperatures; this method updates those scalars immediately before
		graph.launch(), avoiding by-value capture without switching back to the
		legacy pure-eager route.
		"""
		if plan is None:
			return
		integ = plan.get("integ")
		fn = getattr(integ, "graph_prepare_replay", None)
		if callable(fn):
			fn(timestep)
		for b in plan.get("backends", ()): 
			force_graph_prepare_replay(b, timestep)

	def _launch_graph(self, graph, kind, plan=None, timestep=None):
		"""Launch a captured CUDA graph.  Timing mode measures the whole replay
		as a single GPU region; kernels inside replay are not decomposed.

		Important: keep replay on the caller/current stream, as in the original
		multi-force graph path.  stochastic DPD graph replay only needs graph-dynamic inputs (DPD seed scalar)
		updated immediately before replay; forcing a private replay stream can race
		against eager sort/topology or post-step work that still uses the default
		stream in mixed cases such as 13_sort_topology_refresh.
		"""
		cat = "cuda_graph_%s_replay" % kind
		# Dynamic graph inputs are enqueued on the same current stream as graph.launch(),
		# so the seed update remains ordered without changing the stream semantics of
		# existing graph-covered cases.
		self._prepare_graph_inputs(plan, timestep)
		if _t.ENABLED:
			_t.gpu_begin(cat)
			try:
				graph.launch()
			finally:
				_t.gpu_end(cat)
		else:
			graph.launch()

	def _ensure_cupy(self):
		if self.cp is not None:
			return self.available
		try:
			import cupy as cp
			self.cp = cp
			self.available = (hasattr(cp.cuda, "Graph") and hasattr(cp.cuda, "Stream")
			                  and hasattr(cp.cuda.Stream, "begin_capture"))
		except Exception:
			self.cp = False
			self.available = False
		return self.available

	# ---- force plan / eligibility --------------------------------------------
	def _force_backend(self, force_obj):
		# Duck-typed graph-force protocol.  Built-in and future custom
		# force wrappers may expose either .data.backend or .data; the graph
		# planner intentionally avoids concrete class checks here.
		return force_graph_backend_from_wrapper(force_obj)

	def _backend_name(self, backend):
		return force_graph_backend_name(backend)

	def _unready_reason(self, backend):
		return force_graph_unready_reason(backend)

	def _nlist_wrappers(self, backends):
		out = []
		seen = set()
		for b in backends:
			nlw = getattr(b, "nlist", None)
			if nlw is None:
				continue
			key = id(nlw)
			if key not in seen:
				seen.add(key)
				out.append(nlw)
		return out

	def _integrator_unready_reason(self, integ):
		fn = getattr(integ, "graph_unready_reason", None)
		if callable(fn):
			try:
				return fn()
			except Exception:
				return "unsupported_integrator"
		return "unsupported_integrator"

	def _integrator_target(self, integ):
		"""Return a capture-signature target value for graph-compatible integrators.

		RawNVT/RawBD return a stable signature target.  Their captured kernels read
		the live thermostat temperature from a one-element device scalar updated by
		graph_prepare_replay(timestep), so variant temperatures do not force graph
		fallback.  RawNVE has no thermostat/temperature target and returns None.
		"""
		fn = getattr(integ, "graph_target_value", None)
		if callable(fn):
			try:
				return float(fn())
			except Exception:
				self._fb("variant_temp"); return None
		if hasattr(integ, "graph_temp_partial") and hasattr(integ, "graph_xi_update"):
			if not hasattr(integ, "temp") or hasattr(integ.temp, "getValue"):
				self._fb("variant_temp"); return None
			try:
				target = float(integ.temp)
			except Exception:
				self._fb("variant_temp"); return None
			if target <= 0.0:
				self._fb("unsupported_integrator"); return None
			return target
		# NVE-like integrator: no target stage.
		return None

	def _integrator_signature(self, integ):
		fn = getattr(integ, "graph_signature", None)
		if callable(fn):
			try:
				return fn()
			except Exception:
				return ("integrator_signature_failed", integ.__class__.__name__, id(integ))
		# Backward-compatible RawNVT signature fallback.
		return (integ.__class__.__name__, id(getattr(integ, "d_xi", None)),
		        id(getattr(integ, "d_partial", None)), int(getattr(integ, "nme", 0)),
		        float(getattr(getattr(integ, "info", None), "dt", 0.0)))

	def _graph_integrator_tail(self, plan):
		integ = plan["integ"]
		if hasattr(integ, "graph_temp_partial") and hasattr(integ, "graph_xi_update"):
			integ.graph_temp_partial()
			integ.graph_xi_update(plan["target"])
		integ.graph_secondstep()

	def _common(self, app):
		"""Return a graph plan dict or None (with fallback reason recorded)."""
		if len(app.integrations) != 1:
			self._fb("multi_integration"); return None
		if len(app.forces) < 1:
			self._fb("no_force"); return None

		integ = app.integrations[0].data
		if not hasattr(integ, "graph_firststep") or not hasattr(integ, "graph_secondstep") or getattr(integ, "fusion", False):
			self._fb("unsupported_integrator"); return None
		ready_fn = getattr(integ, "graph_ready", None)
		if not callable(ready_fn) or not ready_fn():
			self._fb(self._integrator_unready_reason(integ)); return None
		target = self._integrator_target(integ)
		if target is None and (hasattr(integ, "graph_temp_partial") or hasattr(integ, "graph_xi_update")):
			# Thermostat-like integrator failed to produce a static target.
			return None

		backends = []
		for f in app.forces:
			b = self._force_backend(f)
			if b is None:
				self._fb("unsupported_force_wrapper"); return None
			backends.append(b)

		# Legacy single-backend fast path: a graph_force_write() backend overwrites
		# force directly, so no captured zero_force_vp is needed.  This is no longer
		# tied to a concrete RawPair class; any protocol backend can advertise write.
		single_write = (len(backends) == 1 and force_graph_can_write(backends[0]))
		force_mode = "write" if single_write else "add"

		protocol_summary = [describe_force_graph_backend(b) for b in backends]

		for b in backends:
			# Graph-force path is force-only; if virial/potential is required this step,
			# keep eager so dump/compute_properties stay correct.
			need_full_fn = getattr(b, "_need_virial_potential", None)
			if callable(need_full_fn) and bool(need_full_fn()):
				self._fb("need_full_pair"); return None
			if force_mode == "write":
				if not force_graph_can_write(b):
					self._fb("force_graph_missing_graph_force_write"); return None
			else:
				if not force_graph_can_add(b):
					self._fb("force_graph_missing_graph_force_add"); return None
			if not force_graph_ready(b):
				self._fb(self._unready_reason(b)); return None

		self.force_graph_protocol_summary = protocol_summary

		return {
			"app": app,
			"backends": tuple(backends),
			"backend_names": tuple(self._backend_name(b) for b in backends),
			"force_graph_protocol_summary": tuple(protocol_summary),
			"force_mode": force_mode,
			"integ": integ,
			"target": target,
			"nlists": tuple(self._nlist_wrappers(backends)),
		}

	# ---- FULL graph (firststep included) ----
	def try_run_full(self, app, timestep):
		if not self.enabled or not self.allow_full or self._disabled:
			return False
		_sp_inc("cuda_graph_full_attempt_count")
		if not self._ensure_cupy() or not self.available:
			self._disabled = True; self._fb("unsupported_cuda_graph"); return False
		plan = self._common(app)
		if plan is None:
			return False
		for nlw in plan["nlists"]:
			nl = getattr(nlw, "data", None)
			if nl is None or not hasattr(nl, "will_do_gpu_check"):
				self._fb("full_nlist_check"); return False
			if nl.will_do_gpu_check(False):
				self._fb("full_nlist_check"); return False
		for nlw in plan["nlists"]:
			nlw.calculate(timestep)
			nl = getattr(nlw, "data", None)
			if getattr(nl, "did_gpu_check_this_call", True):
				self._fb("full_nlist_check"); return False
		sig = self._full_signature(plan)
		target = plan["target"]
		if (self.full_graph is None) or (self.full_sig != sig) or (self.full_target != target):
			if not self._capture_full(plan, sig, timestep):
				return False
			self.full_replay += 1
			self.full_step_count += 1
			prof.count("cuda_graph_full_replay_count")
			_sp_inc("cuda_graph_full_replay_count")
			return True
		self._launch_graph(self.full_graph, "full", plan, timestep)
		self.full_replay += 1
		self.full_step_count += 1
		prof.count("cuda_graph_full_replay_count")
		_sp_inc("cuda_graph_full_replay_count")
		return True

	# ---- TAIL graph (eager firststep + eager nlist already done) --------------
	def prepare_tail_nlists(self, app, timestep):
		"""Run authoritative nlist.calculate() calls before tail graph.  Stores the
		validated plan so try_run_tail() can avoid rebuilding it with different state."""
		self._tail_plan = None
		if not self.enabled or self._disabled:
			return False
		if not self._ensure_cupy() or not self.available:
			self._disabled = True; self._fb("unsupported_cuda_graph"); return False
		plan = self._common(app)
		if plan is None:
			return False
		for nlw in plan["nlists"]:
			nlw.calculate(timestep)
		self._tail_plan = plan
		return True

	def try_run_tail(self, app, timestep):
		if not self.enabled or self._disabled:
			return False
		_sp_inc("cuda_graph_tail_attempt_count")
		plan = self._tail_plan
		self._tail_plan = None
		if plan is None:
			plan = self._common(app)
			if plan is None:
				return False
		for nlw in plan["nlists"]:
			nl = getattr(nlw, "data", None)
			if nl is None or getattr(nl, "did_gpu_check_this_call", True):
				self._fb("nlist_check"); return False
		sig = self._tail_signature(plan)
		target = plan["target"]
		if (self.tail_graph is None) or (self.tail_sig != sig) or (self.tail_target != target):
			if not self._capture_tail(plan, sig, timestep):
				return False
			self.tail_replay += 1
			self.tail_step_count += 1
			prof.count("cuda_graph_tail_replay_count")
			_sp_inc("cuda_graph_tail_replay_count")
			return True
		self._launch_graph(self.tail_graph, "tail", plan, timestep)
		self.tail_replay += 1
		self.tail_step_count += 1
		prof.count("cuda_graph_tail_replay_count")
		_sp_inc("cuda_graph_tail_replay_count")
		return True

	# ---- signatures -----------------------------------------------------------
	def _force_signature(self, plan):
		parts = []
		for b in plan["backends"]:
			parts.append((self._backend_name(b), force_graph_signature(b)))
		return tuple(parts)

	def _full_signature(self, plan):
		return ("full", plan["force_mode"], plan["backend_names"], self._force_signature(plan),
		        self._integrator_signature(plan["integ"]),
		        id(getattr(plan["app"].info.cu, "d_force", None)),
		        id(getattr(plan["app"].info.cu, "d_virial_potential", None)))

	def _tail_signature(self, plan):
		return ("tail", plan["force_mode"], plan["backend_names"], self._force_signature(plan),
		        self._integrator_signature(plan["integ"]),
		        id(getattr(plan["app"].info.cu, "d_force", None)),
		        id(getattr(plan["app"].info.cu, "d_virial_potential", None)))

	# ---- capture helpers ------------------------------------------------------
	def _new_stream(self):
		if self.stream is None:
			self.stream = self.cp.cuda.Stream(non_blocking=True)
		return self.stream

	def _graph_zero_force(self, app):
		cu = getattr(app.info, "cu", None)
		if cu is None:
			raise RuntimeError("CUDA Graph additive force segment requires CuInfo arrays")
		cu.d_force.fill(0.0)
		cu.d_virial_potential.fill(0.0)

	def _graph_force_segment(self, plan):
		if plan["force_mode"] == "write":
			plan["backends"][0].graph_force_write()
		else:
			self._graph_zero_force(plan["app"])
			for b in plan["backends"]:
				b.graph_force_add()

	def _capture_full(self, plan, sig, timestep):
		s = self._new_stream()
		integ = plan["integ"]
		target = plan["target"]
		try:
			with s:
				s.begin_capture()
				integ.graph_firststep()
				self._graph_force_segment(plan)
				self._graph_integrator_tail(plan)
				self.full_graph = s.end_capture()
			self.full_sig = sig
			self.full_target = target
			self.full_capture += 1
			prof.count("cuda_graph_full_capture_count")
			_sp_inc("cuda_graph_full_capture_count")
			self._launch_graph(self.full_graph, "full", plan, timestep)  # capture does not execute -> run this step
			return True
		except Exception as e:
			self._disabled = True
			self._fb("capture_failed")
			if self.gprofile:
				print("info : [cuda_graph] full capture failed -> eager: %s" % repr(e))
			try: s.end_capture()
			except Exception: pass
			self.full_graph = None
			return False

	def _capture_tail(self, plan, sig, timestep):
		s = self._new_stream()
		integ = plan["integ"]
		target = plan["target"]
		try:
			with s:
				s.begin_capture()
				self._graph_force_segment(plan)
				self._graph_integrator_tail(plan)
				self.tail_graph = s.end_capture()
			self.tail_sig = sig
			self.tail_target = target
			self.tail_capture += 1
			prof.count("cuda_graph_tail_capture_count")
			_sp_inc("cuda_graph_tail_capture_count")
			self._launch_graph(self.tail_graph, "tail", plan, timestep)
			return True
		except Exception as e:
			self._disabled = True
			self._fb("capture_failed")
			if self.gprofile:
				print("info : [cuda_graph] tail capture failed -> eager: %s" % repr(e))
			try: s.end_capture()
			except Exception: pass
			self.tail_graph = None
			return False

	def diag_snapshot(self):
		return self.diagnostics()

	def _coverage(self):
		total = int(self.full_step_count) + int(self.tail_step_count) + int(self.eager_step_count)
		graph = int(self.full_step_count) + int(self.tail_step_count)
		return (float(graph) / float(total)) if total > 0 else None

	def diagnostics(self):
		cov = self._coverage()
		loop_kind = getattr(self, "application_loop_kind", getattr(self, "application_loop", "unknown"))
		loop_label = getattr(self, "application_loop_label", loop_kind)
		loop_count = int(getattr(self, "application_loop_count", 0) or 0)
		loop_eager_graph_count = int(getattr(self, "application_loop_eager_graph_count", 0) or 0)
		loop_eager_debug_count = int(getattr(self, "application_loop_eager_debug_count", 0) or 0)
		loop_graph_aware = bool(loop_kind != "eager_debug" and loop_kind not in ("unknown", "not_started"))
		out = {
			"pygamd_execution_policy": self.execution_policy,
			"pygamd_execution_policy_explicit": bool(getattr(self.policy, "explicit", False)),
			"pygamd_execution_route_tier": getattr(self.policy, "route_tier", "unknown"),
			"pygamd_execution_route_label": getattr(self.policy, "route_label", self.execution_policy),
			"pygamd_execution_is_primary_route": bool(getattr(self.policy, "is_primary_route", False)),
			"pygamd_execution_is_legacy_debug_route": bool(getattr(self.policy, "is_legacy_debug_route", False)),
			"force_graph_protocol_version": self.force_graph_protocol_version,
			"force_graph_protocol_summary": list(self.force_graph_protocol_summary),
			"force_graph_protocol_backend_names": [x.get("backend_name") for x in self.force_graph_protocol_summary],
			"pygamd_application_loop_kind": loop_kind,
			"pygamd_application_loop_label": loop_label,
			"pygamd_application_loop_count": loop_count,
			"pygamd_application_loop_eager_graph_count": loop_eager_graph_count,
			"pygamd_application_loop_eager_debug_count": loop_eager_debug_count,
			# Legacy aliases retained for report compatibility.
			"pygamd_application_loop": loop_kind,
			"pygamd_application_loop_graph_aware": loop_graph_aware,
			"pygamd_application_loop_is_primary": loop_kind == "eager_graph_primary",
			"pygamd_application_loop_is_eager_debug": loop_kind == "eager_debug",
			"cuda_graph_enabled": bool(self.enabled),
			"cuda_graph_mode": self.mode,
			"cuda_graph_strict": bool(self.strict),
			"cuda_graph_strict_min_coverage": float(self.strict_min_coverage),
			"cuda_graph_strict_max_eager_steps": self.strict_max_eager_steps,
			"cuda_graph_diag_enabled": bool(self.gdiag),
			"cuda_graph_diag_env_value": os.environ.get("PYGAMD_CUDA_GRAPH_DIAG", ""),
			"cuda_graph_full_capture_count": int(self.full_capture),
			"cuda_graph_full_replay_count": int(self.full_replay),
			"cuda_graph_tail_capture_count": int(self.tail_capture),
			"cuda_graph_tail_replay_count": int(self.tail_replay),
			"cuda_graph_full_step_count": int(self.full_step_count),
			"cuda_graph_tail_step_count": int(self.tail_step_count),
			"cuda_graph_eager_step_count": int(self.eager_step_count),
			"cuda_graph_graph_step_count": int(self.full_step_count) + int(self.tail_step_count),
			"cuda_graph_total_policy_step_count": int(self.full_step_count) + int(self.tail_step_count) + int(self.eager_step_count),
			"cuda_graph_coverage_pct": (cov * 100.0 if cov is not None else None),
			"cuda_graph_fallback_count": int(self.fallback_count),
			"cuda_graph_fallback_reasons": dict(self.reasons),
			"cuda_graph_disabled_internal": bool(self._disabled),
		}
		return out

	def _check_strict(self):
		if not self.strict:
			return
		total = int(self.full_step_count) + int(self.tail_step_count) + int(self.eager_step_count)
		if total <= 0:
			return
		coverage = self._coverage()
		if self.strict_max_eager_steps is not None and int(self.eager_step_count) > int(self.strict_max_eager_steps):
			raise RuntimeError(
				"PYGAMD graph_strict failed: eager fallback steps %d exceed max %d; reasons=%s"
				% (self.eager_step_count, self.strict_max_eager_steps, dict(self.reasons)))
		if coverage is not None and coverage < float(self.strict_min_coverage):
			raise RuntimeError(
				"PYGAMD graph_strict failed: graph coverage %.2f%% < %.2f%%; graph_steps=%d eager_steps=%d reasons=%s"
				% (coverage * 100.0, float(self.strict_min_coverage) * 100.0,
				   int(self.full_step_count) + int(self.tail_step_count),
				   int(self.eager_step_count), dict(self.reasons)))

	def finalize(self):
		if self.gprofile or self.strict:
			print("info : [cuda_graph] policy=%s mode=%s full(cap=%d,replay=%d,steps=%d) tail(cap=%d,replay=%d,steps=%d) eager_steps=%d coverage=%s fallback=%d reasons=%s"
			      % (self.execution_policy, self.mode, self.full_capture, self.full_replay, self.full_step_count,
			         self.tail_capture, self.tail_replay, self.tail_step_count, self.eager_step_count,
			         ("%.2f%%" % (self._coverage() * 100.0) if self._coverage() is not None else "n/a"),
			         self.fallback_count, dict(self.reasons)))
		self._check_strict()
