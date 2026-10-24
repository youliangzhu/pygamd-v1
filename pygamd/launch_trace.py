'''
Launch trace / step-loop attribution (opt-in, ZERO overhead when disabled).

Enable with PYGAMD_LAUNCH_TRACE=1.  This is a SEPARATE switch from the
engine timing layer (pygamd.engine_timing / PYGAMD_ENGINE_TIMING): launch
tracing only COUNTS events (kernel launches, per-run actions, steps) with plain
dict increments -- it performs NO device synchronize, NO print, NO file write, and
makes no allocation in the hot path.  When unset/0/false/off every entry point is
a single bool check + return, and call sites are additionally guarded by
`if _lt.LAUNCH:` so the default run path is unchanged.

Why separate from timing: launch counts are useful WITHOUT the per-region
synchronize that timing adds, so we can attribute "how many kernels per step"
without the timing perturbation.  Combine with PYGAMD_ENGINE_TIMING to also get
per-category / per-stage GPU ms (the launch collector merges both).

Imports only os at module load (no cupy, no numba), so `import pygamd` stays
CUDA-driver-free.  This module measures NOTHING in ms; ms come from engine_timing.
'''

import os


def _truthy(v):
    return v not in ("", "0", "false", "False", "no", "off", "NO", "OFF", None)


# --- env + sentinel force-enable (mirrors the runtime timing mechanism) ----------
# LAUNCH is normally driven by PYGAMD_LAUNCH_TRACE, read at IMPORT time.
# A driver may import pygamd before the env is set, or a subprocess may not
# inherit it, so the launch trace runner also drops "force launch trace" SENTINEL
# files on disk; if any is present launch tracing is enabled.  force_enable()
# re-evaluates env+sentinels and can flip LAUNCH ON later (NEVER OFF).
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))   # .../pygamd
_ROOT = os.path.dirname(_PKG_DIR)                         # project root


def _launch_flag_paths():
    paths = [
        os.path.join(_PKG_DIR, "_pygamd_force_launch_trace.flag"),
        os.path.join(_ROOT, ".pygamd_force_launch_trace.flag"),
        "/tmp/pygamd_force_launch_trace.flag",
    ]
    try:
        paths.append("/tmp/pygamd_force_launch_trace_%d.flag" % os.getuid())
    except Exception:
        pass
    return paths


def _launch_flag_present():
    for p in _launch_flag_paths():
        try:
            if os.path.exists(p):
                return True
        except Exception:
            pass
    return False


def _launch_env_truthy():
    return _truthy(os.environ.get("PYGAMD_LAUNCH_TRACE", ""))


LAUNCH_TRACE_ENV_VALUE = os.environ.get("PYGAMD_LAUNCH_TRACE", "")
LAUNCH = _launch_env_truthy() or _launch_flag_present()


def force_enable():
    '''Re-evaluate LAUNCH from PYGAMD_LAUNCH_TRACE + sentinel files.
    Idempotent; NEVER turns tracing OFF.  Lets the runner enable launch tracing
    even if this module was imported before the env / sentinel existed.'''
    global LAUNCH, LAUNCH_TRACE_ENV_VALUE
    LAUNCH_TRACE_ENV_VALUE = os.environ.get("PYGAMD_LAUNCH_TRACE", "")
    if not LAUNCH and (_launch_env_truthy() or _launch_flag_present()):
        LAUNCH = True
    return LAUNCH

# kernel-launch categories (one increment per real RawKernel launch).
LAUNCH_CATEGORIES = (
    "nlist_check_kernel", "clist_build_kernel", "nlist_build_kernel",
    "pair_force", "dpd_force", "bond_force",
    "angle_force", "dihedral_force",
    # NVT decomposition (source-accurate kernel stages, see integrations/nvt.py)
    "integrator_nvt_firststep", "integrator_nvt_ke_reduce",
    "integrator_nvt_xi_secondstep", "integrator_nvt_xi_update",
    "integrator_nvt_secondstep",
    "integrator_nve", "integrator_bd",
    # host-side "launch-like" calls (not GPU kernels, counted for completeness)
    "sort_launches_or_host_calls", "topology_refresh_calls",
)
# NVT stage launches summed into the derived integrator_nvt_launches total.
NVT_STAGE_LAUNCHES = (
    "integrator_nvt_firststep", "integrator_nvt_ke_reduce",
    "integrator_nvt_xi_secondstep", "integrator_nvt_xi_update",
    "integrator_nvt_secondstep",
)
# per-run action counters (things the step loop does, not all GPU kernels).
ACTION_NAMES = (
    "register_calls", "force_compute_calls", "dump_calls",
    "integrator_firststep_calls", "integrator_secondstep_calls",
)

_launches = {}     # name -> count
_actions = {}      # name -> count
_step = 0


def launch(name):
    '''Count one real kernel launch for `name` (no-op when tracing is off).'''
    if not LAUNCH:
        return
    _launches[name] = _launches.get(name, 0) + 1


def action(name):
    '''Count one per-run action (register/force-dispatch/dump/integrator call).'''
    if not LAUNCH:
        return
    _actions[name] = _actions.get(name, 0) + 1


def tick_step():
    '''Count one Python step-loop iteration (no-op when tracing is off).'''
    if not LAUNCH:
        return
    global _step
    _step += 1


def reset():
    '''Clear counters/step (call after warmup, before the timed segment).'''
    _launches.clear()
    _actions.clear()
    global _step
    _step = 0


def snapshot():
    '''Return a flat dict of <category>_launches + per-run actions + totals.
    Returns {"launch_trace_enabled": False} when tracing is off.'''
    if not LAUNCH:
        return {"launch_trace_enabled": False}
    out = {"launch_trace_enabled": True,
           "launch_force_flag_active": _launch_flag_present(),
           "launch_trace_env_value": LAUNCH_TRACE_ENV_VALUE}
    total = 0
    for cat in LAUNCH_CATEGORIES:
        c = _launches.get(cat, 0)
        out[cat + "_launches"] = c
        total += c
    # derived NVT total launches = sum of its stages (display; already in total)
    out["integrator_nvt_launches"] = sum(_launches.get(s, 0) for s in NVT_STAGE_LAUNCHES)
    out["launch_count_total"] = total
    out["step_count"] = _step
    out["launch_count_per_step"] = round(total / _step, 4) if _step else None
    for a in ACTION_NAMES:
        out[a] = _actions.get(a, 0)
    return out
