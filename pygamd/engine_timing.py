'''
Engine timing decomposition (opt-in, ZERO overhead when disabled).

Enable with PYGAMD_ENGINE_TIMING=1.  When unset/0/false/off, every entry point is
a single module-global bool check + return, and call sites are additionally
guarded by `if _t.ENABLED:` so the default run path produces NO timing work at
all (no CUDA events, no synchronize, no allocation).

This is a BOTTLENECK-ATTRIBUTION layer, not an optimizer.  It answers "where does
the time go" per big performance case:

  * GPU kernel regions  -> measured with CuPy CUDA Events (gpu_begin/gpu_end).
    get_elapsed_time gives GPU-only elapsed ms.  gpu_end SYNCHRONIZES that region
    so the elapsed can be read -- this is a deliberate measurement perturbation
    that ONLY happens in timing mode (documented in the report).
  * Host copy / sync points (cp.asnumpy, reallocation, host-key sort) -> measured
    with Python wall time (host_begin/host_end).  These are real synchronization
    points in the normal code, so wall time there is fair.

Aggregation only: counts/ms accumulate per category; nothing is printed per step
(unless PYGAMD_ENGINE_TIMING_VERBOSE=1, which prints a one-line summary every
PYGAMD_ENGINE_TIMING_EVERY steps -- default 100, never on by default).

Imports only os/time at module load (no cupy, no numba), so `import pygamd` stays
CUDA-driver-free; cupy is imported lazily on first GPU measurement.
'''

import os
import time


def _truthy(v):
    return v not in ("", "0", "false", "False", "no", "off", "NO", "OFF", None)


# --- env + sentinel force-enable (runtime timing mechanism) ------------------
# ENABLED is normally driven by PYGAMD_ENGINE_TIMING.  But the env var is read at
# IMPORT time, and a driver may import pygamd before the env is set, or a
# subprocess may not inherit it.  As a belt-and-suspenders the engine metrics runners
# also drop "force timing" SENTINEL files on disk; if any is present timing is
# enabled.  force_enable() re-evaluates env+sentinels and can flip ENABLED ON
# later (it NEVER turns OFF), so reset()/snapshot() can fix an early import.
_PKG_DIR = os.path.dirname(os.path.abspath(__file__))   # .../pygamd
_ROOT = os.path.dirname(_PKG_DIR)                         # project root


def _timing_flag_paths():
    paths = [
        os.path.join(_PKG_DIR, "_pygamd_force_timing.flag"),
        os.path.join(_ROOT, ".pygamd_force_timing.flag"),
        "/tmp/pygamd_force_timing.flag",
    ]
    try:
        paths.append("/tmp/pygamd_force_timing_%d.flag" % os.getuid())
    except Exception:
        pass
    return paths


def _timing_flag_present():
    for p in _timing_flag_paths():
        try:
            if os.path.exists(p):
                return True
        except Exception:
            pass
    return False


def _timing_env_truthy():
    return _truthy(os.environ.get("PYGAMD_ENGINE_TIMING", ""))


TIMING_ENV_VALUE = os.environ.get("PYGAMD_ENGINE_TIMING", "")
ENABLED = _timing_env_truthy() or _timing_flag_present()
VERBOSE = _truthy(os.environ.get("PYGAMD_ENGINE_TIMING_VERBOSE", ""))
try:
    EVERY = int(os.environ.get("PYGAMD_ENGINE_TIMING_EVERY", "100"))
except (TypeError, ValueError):
    EVERY = 100
if EVERY < 1:
    EVERY = 100


def force_enable():
    '''Re-evaluate ENABLED from PYGAMD_ENGINE_TIMING + sentinel files.  Idempotent;
    NEVER turns timing OFF.  Lets the runner enable timing even if this module was
    imported before the env / sentinel existed.  Returns the (possibly updated)
    ENABLED.'''
    global ENABLED, TIMING_ENV_VALUE
    TIMING_ENV_VALUE = os.environ.get("PYGAMD_ENGINE_TIMING", "")
    if not ENABLED and (_timing_env_truthy() or _timing_flag_present()):
        ENABLED = True
    return ENABLED

# GPU kernel leaf categories (measured with CUDA events; summed into the total).
# launch trace: integrator_nvt is now decomposed into source-accurate kernel stages
# (NVT_STAGE_CATEGORIES); integrator_nvt_ms is DERIVED as the sum of those stages
# so the timing report keeps working without double-counting.
GPU_CATEGORIES = (
    "nlist_check_kernel", "clist_build_kernel", "nlist_build_kernel",
    "pair_force", "dpd_force", "bond_force",
    "angle_force", "dihedral_force",
    "integrator_nve", "integrator_bd",
    # CUDA Graph diagnostics: graph replay is timed as one whole-graph GPU region.
    # The kernels *inside* the graph are not decomposed here.
    "cuda_graph_full_replay", "cuda_graph_tail_replay",
)
NVT_STAGE_CATEGORIES = (
    "integrator_nvt_firststep", "integrator_nvt_ke_reduce",
    "integrator_nvt_xi_secondstep", "integrator_nvt_xi_update",
    "integrator_nvt_secondstep",
)
# Host copy / sync categories (measured with Python wall time).
HOST_CATEGORIES = (
    "nlist_displacement_copy", "clist_status_copy", "nlist_width_copy",
    "nlist_fused_status_copy", "nlist_realloc", "clist_realloc",
    "sort", "topology_refresh",
)
# Categories we do NOT instrument in this patch -> reported as not_available.
# CUDA Graph diagnostics instruments whole graph replay, but not graph capture nor kernels
# inside a replayed graph.
NOT_INSTRUMENTED = ("cuda_graph_capture", "cuda_graph_internal_kernels")

PERTURBATION_NOTE = (
    "timing mode records CUDA events around each GPU region and SYNCHRONIZES at "
    "region end to read GPU-only elapsed ms; host copies use wall time.  The "
    "per-region synchronize is a measurement perturbation that inflates absolute "
    "step time -- use the PROPORTIONS here, not the absolute TPS, from a timing run."
)

_cp = None                 # lazy cupy handle (None=unknown, False=absent)
_gpu = {}                  # name -> [calls, total_ms]
_host = {}                 # name -> [calls, total_ms]
_events = {}               # name -> (start_event, stop_event)  (reused pool)
_gpu_open = {}             # name -> (start_event, stop_event)  (currently open)
_host_open = {}            # name -> perf_counter t0
_step = 0


def _get_cp():
    global _cp
    if _cp is False:
        return None
    if _cp is None:
        try:
            import cupy as cp
            _cp = cp
        except Exception:
            _cp = False
            return None
    return _cp


# --- GPU region timing (CUDA events) -----------------------------------------
def gpu_begin(name):
    '''Record a start CUDA event for `name` on the current stream (no sync).'''
    if not ENABLED:
        return
    cp = _get_cp()
    if cp is None:
        return
    ev = _events.get(name)
    if ev is None:
        ev = (cp.cuda.Event(), cp.cuda.Event())
        _events[name] = ev
    ev[0].record()
    _gpu_open[name] = ev


def gpu_end(name):
    '''Record the stop event, SYNCHRONIZE it (timing-mode perturbation), and
    accumulate GPU-only elapsed ms for `name`.'''
    if not ENABLED:
        return
    ev = _gpu_open.pop(name, None)
    if ev is None:
        return
    cp = _get_cp()
    if cp is None:
        return
    ev[1].record()
    ev[1].synchronize()
    dt = float(cp.cuda.get_elapsed_time(ev[0], ev[1]))   # ms (GPU-only)
    e = _gpu.get(name)
    if e is None:
        _gpu[name] = [1, dt]
    else:
        e[0] += 1
        e[1] += dt


# --- host copy / sync timing (wall time) -------------------------------------
def host_begin(name):
    if not ENABLED:
        return
    _host_open[name] = time.perf_counter()


def host_end(name):
    if not ENABLED:
        return
    t0 = _host_open.pop(name, None)
    if t0 is None:
        return
    dt = (time.perf_counter() - t0) * 1000.0             # ms (wall)
    e = _host.get(name)
    if e is None:
        _host[name] = [1, dt]
    else:
        e[0] += 1
        e[1] += dt


# --- per-step hook (verbose only) --------------------------------------------
def tick_step():
    if not ENABLED:
        return
    global _step
    _step += 1
    if VERBOSE and (_step % EVERY == 0):
        top = max(list(_gpu.items()) + list(_host.items()),
                  key=lambda kv: kv[1][1], default=None)
        if top is not None:
            print("info : [engine timing] step=%d top=%s %.3f ms"
                  % (_step, top[0], top[1][1]))


def reset():
    '''Clear accumulated timings/step (keep the event pool).  Call after warmup,
    before the timed segment.'''
    _gpu.clear()
    _host.clear()
    _gpu_open.clear()
    _host_open.clear()
    global _step
    _step = 0


def snapshot(wall_timed_region_ms=None):
    '''Return a flat dict of <category>_ms / <category>_calls plus totals.

    Instrumented-but-unused categories report 0.0 (the case did not hit them);
    NOT_INSTRUMENTED categories report "not_available" with a reason.  Returns a
    minimal {"timing_enabled": False} when timing is off.'''
    if not ENABLED:
        return {"timing_enabled": False}
    out = {"timing_enabled": True, "timing_perturbation_note": PERTURBATION_NOTE,
           "timing_force_flag_active": _timing_flag_present(),
           "timing_env_value": TIMING_ENV_VALUE}
    total = 0.0
    gpu_total = 0.0
    gpu_calls = 0
    for cat in GPU_CATEGORIES + NVT_STAGE_CATEGORIES:
        e = _gpu.get(cat)
        out[cat + "_ms"] = round(e[1], 4) if e else 0.0
        out[cat + "_calls"] = e[0] if e else 0
        if e:
            total += e[1]
            gpu_total += e[1]
            gpu_calls += e[0]
    # Derived coarse integrator_nvt (sum of stages) for the timing report; NOT added
    # to the total again (its stages are already counted above).
    nvt_ms = sum(_gpu[s][1] for s in NVT_STAGE_CATEGORIES if s in _gpu)
    nvt_calls = sum(_gpu[s][0] for s in NVT_STAGE_CATEGORIES if s in _gpu)
    out["integrator_nvt_ms"] = round(nvt_ms, 4)
    out["integrator_nvt_calls"] = nvt_calls
    host_total = 0.0
    for cat in HOST_CATEGORIES:
        e = _host.get(cat)
        out[cat + "_ms"] = round(e[1], 4) if e else 0.0
        out[cat + "_calls"] = e[0] if e else 0
        if e:
            total += e[1]
            host_total += e[1]
    for cat in NOT_INSTRUMENTED:
        out[cat + "_ms"] = "not_available"
    # CUDA Graph diagnostics: whole graph replay is now timed around graph.launch().
    # This splits graph replay from unattributed_ms while still NOT decomposing
    # the kernels inside the graph.
    full_g = _gpu.get("cuda_graph_full_replay")
    tail_g = _gpu.get("cuda_graph_tail_replay")
    graph_ms = (full_g[1] if full_g else 0.0) + (tail_g[1] if tail_g else 0.0)
    graph_calls = (full_g[0] if full_g else 0) + (tail_g[0] if tail_g else 0)
    out["cuda_graph_replay_ms"] = round(graph_ms, 4)
    out["cuda_graph_replay_calls"] = graph_calls
    out["cuda_graph_replay_timing_mode"] = "whole_graph_replay_only"
    out["gpu_profiled_ms"] = round(gpu_total, 4)
    out["host_profiled_ms"] = round(host_total, 4)
    out["gpu_region_count_total"] = gpu_calls
    out["cuda_graph_timing_reason"] = (
        "CUDA Graph replay is timed as a whole graph.launch() region when "
        "PYGAMD_ENGINE_TIMING is active; kernels inside the replayed graph are "
        "not decomposed by this timing layer.")
    out["other_profiled_ms"] = 0.0          # catch-all; nothing falls here yet
    out["total_profiled_ms"] = round(total, 4)
    if isinstance(wall_timed_region_ms, (int, float)):
        out["wall_timed_region_ms"] = round(float(wall_timed_region_ms), 4)
        out["unattributed_ms"] = round(max(float(wall_timed_region_ms) - total, 0.0), 4)
    else:
        out["wall_timed_region_ms"] = "not_available"
        out["unattributed_ms"] = "not_available"
    return out
