'''
Engine metrics counters (opt-in, near-zero overhead when disabled).

Enable with the environment variable PYGAMD_ENGINE_METRICS=1.  When disabled,
every entry point is a single module-global bool check + return, so the hot path
is unaffected (this module imports only os/json/time; no cupy, no numba, so
`import pygamd` stays CUDA-driver-free).

This is a LIGHTWEIGHT COUNTER layer (it does NOT add device synchronizations the
way pygamd.profiler does).  It counts events per run so engine metrics can locate the
bottleneck (Python step loop / kernel launches / nlist / clist / force kernels /
sort host-copies) without changing physics or default behavior.

Typical usage inside instrumented modules:
    from pygamd import engine_metrics as _sp
    _sp.inc("nlist_check_count")

The engine metrics runner controls run metadata + reset:
    import pygamd.engine_metrics as sp
    sp.reset()                      # after warmup, before the timed segment
    ... app.run(timed) ...
    sp.record_graph(app)            # pull CUDA-graph capture/replay counts
    sp.dump_json("engine_metrics.json",
                 extra=dict(case=..., size_profile=..., engine="cupy",
                            warmup=..., timed=..., total_wall_s=..., internal_tps=...))
'''

import os
import json
import time

from pygamd import engine_timing as _timing   # engine timing opt-in timing layer
from pygamd import launch_trace as _ltrace  # launch trace opt-in launch trace

# Read once at import.  The engine metrics runner sets the env BEFORE importing pygamd.
ENABLED = os.environ.get("PYGAMD_ENGINE_METRICS", "") not in (
    "", "0", "false", "False", "no", "off", "NO", "OFF")

_counters = {}
_times = {}
_notes = []
_meta = {}

# canonical counter names engine metrics reports on (kept here for the collector / docs)
COUNTER_NAMES = (
    "step_loop_iters",
    "nlist_check_count", "nlist_rebuild_count", "clist_build_count",
    "nlist_width_realloc_count", "force_clear_count",
    "pair_compute_count", "pair_charge_compute_count", "dpd_compute_count", "pppm_compute_count",
    "bond_compute_count", "angle_compute_count", "dihedral_compute_count",
    "integration_firststep_count", "integration_secondstep_count",
    "sort_count", "host_sync_count",
    # neighbor-list diagnostics: nlist rebuild-reason breakdown (sums to nlist_rebuild_count;
    # width_realloc counts EXTRA build iterations, so it is additive on top).
    "nlist_rebuild_due_to_periodic_check",
    "nlist_rebuild_due_to_displacement",
    "nlist_rebuild_due_to_width_realloc",
    "nlist_rebuild_due_to_forced",
    # neighbor-list diagnostics: host-sync reason breakdown (sums to host_sync_count).
    "host_sync_reason_displacement_check",
    "host_sync_reason_width_check",
    "host_sync_reason_clist",
    "host_sync_reason_sort",
    "host_sync_reason_unknown",
    # fused-status: opt-in fused build-status read (PYGAMD_NLIST_FUSED_STATUS).
    # When enabled, the clist + neighbor-build status reads share ONE host copy,
    # so host_sync_reason_fused_status replaces clist + width_check per rebuild.
    "host_sync_reason_fused_status",
    "fused_status_copy_count",
    # CUDA Graph diagnostics: CUDA Graph path diagnostics (per timed segment).
    "cuda_graph_full_attempt_count", "cuda_graph_tail_attempt_count",
    "cuda_graph_full_capture_count", "cuda_graph_tail_capture_count",
    "cuda_graph_full_replay_count", "cuda_graph_tail_replay_count",
    "cuda_graph_full_step_count", "cuda_graph_tail_step_count",
    "cuda_graph_eager_step_count", "cuda_graph_fallback_count",
)

# neighbor-list diagnostics/2b meta fields surfaced at top level of the snapshot (may be null
# when the relevant code path did not run / profiling did not capture them).
_META_FIELDS = (
    "nlist_skin_threshold",            # rbuff/2 displacement threshold (set once)
    "nlist_max_displacement_seen",     # max per-check displacement observed (opt-in kernel)
    "nlist_check_period_effective",    # the period actually used after validation
    "nlist_last_check_step",           # timestep/call of the most recent GPU check
    "nlist_last_rebuild_step",         # timestep/call of the most recent rebuild
    # neighbor-list skin policy: opt-in rbuff/skin policy (PYGAMD_NLIST_RBUFF_SCALE).
    "nlist_rbuff_scale_effective",     # rbuff/rcut actually used (default 0.1)
    "nlist_rbuff_scale_source",        # "default" | "env" | "fallback"
    # fused-status: opt-in fused build-status read (PYGAMD_NLIST_FUSED_STATUS).
    "nlist_fused_status_enabled",      # bool: fused clist+width status copy active
    # RawNVT fusion: opt-in RawNVT finalize_xi+secondstep fusion (PYGAMD_RAWNVT_FUSION).
    "rawnvt_fusion_enabled",           # bool: NVT fused path active (4->3 launches)
    "rawnvt_fusion_mode",              # "off" | "safe" | "aggressive"
    "rawnvt_fusion_force_flag_active", # bool: a RawNVT-fusion sentinel file was present
    # CUDA Graph diagnostics: Application-loop and CUDA Graph object-level diagnostics.
    # Canonical application-loop fields are kind/label/count.  Legacy aliases remain below.
    "pygamd_application_loop_kind",
    "pygamd_application_loop_label",
    "pygamd_application_loop_count",
    "pygamd_application_loop_eager_graph_count",
    "pygamd_application_loop_eager_debug_count",
    "pygamd_application_loop",
    "pygamd_application_loop_graph_aware",
    "pygamd_application_loop_is_primary",
    "pygamd_application_loop_is_eager_debug",
    "force_graph_protocol_version",
    "force_graph_protocol_summary",
    "force_graph_protocol_backend_names",
    "cuda_graph_diag_env_value",
    "cuda_graph_diag_enabled",
    "cuda_graph_available",
    "cuda_graph_disabled_internal",
    "cuda_graph_mode",
    "cuda_graph_allow_full",
    "cuda_graph_full_capture_total",
    "cuda_graph_full_replay_total",
    "cuda_graph_tail_capture_total",
    "cuda_graph_tail_replay_total",
    "cuda_graph_fallback_total",
    "cuda_graph_fallback_reasons_total",
)

# counters summed into the conservative kernel_launch_estimate (per-run)
_LAUNCH_COUNTERS = (
    "pair_compute_count", "pair_charge_compute_count", "dpd_compute_count", "pppm_compute_count",
    "bond_compute_count", "angle_compute_count", "dihedral_compute_count",
    "integration_firststep_count", "integration_secondstep_count",
    "clist_build_count", "nlist_rebuild_count", "sort_count",
)


def inc(name, n=1):
    '''Increment a counter (no-op when profiling is disabled).'''
    if not ENABLED:
        return
    _counters[name] = _counters.get(name, 0) + n


def add_time(name, seconds):
    '''Accumulate a timing bucket in seconds (no-op when disabled).'''
    if not ENABLED:
        return
    _times[name] = _times.get(name, 0.0) + float(seconds)


def set_meta(key, value):
    '''Record a scalar metadata value (no-op when disabled).'''
    if not ENABLED:
        return
    _meta[key] = value


def set_meta_max(key, value):
    '''Record the running MAX of a scalar (no-op when disabled).  Used for
    nlist_max_displacement_seen so we keep the largest displacement over the run
    without an extra host sync.'''
    if not ENABLED:
        return
    cur = _meta.get(key)
    if cur is None or value > cur:
        _meta[key] = value


def note(msg):
    '''Record a one-line diagnostic note (deduplicated; no-op when disabled).'''
    if not ENABLED:
        return
    if msg not in _notes:
        _notes.append(msg)


def counter(name):
    return int(_counters.get(name, 0))


def reset():
    '''Clear all counters/times (meta/notes kept).  Call between warmup and the
    timed segment so counts reflect the timed run only.'''
    # engine timing/2e: re-evaluate the opt-in timing / launch-trace switches from
    # env + sentinel files BEFORE the timed segment, so an early `import pygamd`
    # (before the env/sentinel existed) does not leave them stuck OFF.
    _timing.force_enable()
    _ltrace.force_enable()
    _counters.clear()
    _times.clear()
    # keep notes + config meta across reset (they describe the configuration,
    # not the run).  Per-run meta fields are cleared so they reflect the timed
    # segment only (config meta like skin_threshold / check_period_effective stay).
    for k in ("nlist_max_displacement_seen", "nlist_last_check_step",
              "nlist_last_rebuild_step"):
        _meta.pop(k, None)
    # engine timing: clear timing accumulators so they reflect the timed segment only
    # (no-op when PYGAMD_ENGINE_TIMING is off).
    _timing.reset()
    # launch trace: clear launch-trace counters (no-op when launch trace is off).
    _ltrace.reset()


def record_graph(app):
    '''Pull CUDA-graph status + capture/replay counts off a dynamics app.'''
    if not ENABLED:
        return
    g = getattr(app, "_graph", None)
    if g is None:
        set_meta("cuda_graph_enabled", False)
        return
    set_meta("cuda_graph_enabled", bool(getattr(g, "enabled", False)))
    set_meta("cuda_graph_mode", getattr(g, "mode", None))
    # Object totals are cumulative over the app lifetime (useful for sanity),
    # while the engine metrics counters below are per timed segment because reset() clears
    # _counters between warmup and timed runs.
    if hasattr(g, "diag_snapshot"):
        try:
            for k, v in g.diag_snapshot().items():
                set_meta(k, v)
        except Exception as e:
            note("cuda_graph diag_snapshot failed: %r" % (e,))
    obj_cap = int(getattr(g, "full_capture", 0)) + int(getattr(g, "tail_capture", 0))
    obj_rep = int(getattr(g, "full_replay", 0)) + int(getattr(g, "tail_replay", 0))
    # Preserve the legacy top-level fields, but prefer per-run counters when the
    # new CUDA Graph instrumentation populated them.
    run_cap = int(_counters.get("cuda_graph_full_capture_count", 0)) + int(_counters.get("cuda_graph_tail_capture_count", 0))
    run_rep = int(_counters.get("cuda_graph_full_replay_count", 0)) + int(_counters.get("cuda_graph_tail_replay_count", 0))
    _counters["cuda_graph_capture_count"] = run_cap if run_cap else obj_cap
    _counters["cuda_graph_replay_count"] = run_rep if run_rep else obj_rep
    _counters["cuda_graph_capture_count_timed"] = run_cap
    _counters["cuda_graph_replay_count_timed"] = run_rep


def snapshot():
    '''Return a dict of all recorded counters/times/notes/meta.'''
    # engine timing/2e: ensure the opt-in switches reflect env + sentinels at dump
    # time too (idempotent; never turns them off).
    _timing.force_enable()
    _ltrace.force_enable()
    snap = {"counters": dict(_counters), "times_s": dict(_times),
            "notes": list(_notes), "meta": dict(_meta)}
    # surface the canonical counters at top level too (null when not seen)
    for nm in COUNTER_NAMES:
        snap[nm] = _counters.get(nm, None)
    snap["cuda_graph_capture_count"] = _counters.get("cuda_graph_capture_count", None)
    snap["cuda_graph_replay_count"] = _counters.get("cuda_graph_replay_count", None)
    snap["cuda_graph_capture_count_timed"] = _counters.get("cuda_graph_capture_count_timed", None)
    snap["cuda_graph_replay_count_timed"] = _counters.get("cuda_graph_replay_count_timed", None)
    snap["cuda_graph_full_capture_count"] = _counters.get("cuda_graph_full_capture_count", None)
    snap["cuda_graph_tail_capture_count"] = _counters.get("cuda_graph_tail_capture_count", None)
    snap["cuda_graph_full_replay_count"] = _counters.get("cuda_graph_full_replay_count", None)
    snap["cuda_graph_tail_replay_count"] = _counters.get("cuda_graph_tail_replay_count", None)
    snap["cuda_graph_full_step_count"] = _counters.get("cuda_graph_full_step_count", None)
    snap["cuda_graph_tail_step_count"] = _counters.get("cuda_graph_tail_step_count", None)
    snap["cuda_graph_eager_step_count"] = _counters.get("cuda_graph_eager_step_count", None)
    snap["cuda_graph_fallback_count"] = _counters.get("cuda_graph_fallback_count", None)
    snap["cuda_graph_fallback_reasons"] = {
        k.replace("cuda_graph_fallback_reason_", "", 1): v
        for k, v in _counters.items() if k.startswith("cuda_graph_fallback_reason_")
    }
    snap["cuda_graph_enabled"] = _meta.get("cuda_graph_enabled", None)
    # CUDA Graph diagnostics/2k graph diagnostic metadata (recorded at run end by record_graph).
    for nm in (
        "pygamd_execution_policy", "pygamd_execution_policy_explicit",
        "pygamd_execution_route_tier", "pygamd_execution_route_label",
        "pygamd_execution_is_primary_route", "pygamd_execution_is_legacy_debug_route",
        "pygamd_application_loop_kind", "pygamd_application_loop_label",
        "pygamd_application_loop_count", "pygamd_application_loop_eager_graph_count",
        "pygamd_application_loop_eager_debug_count",
        "pygamd_application_loop", "pygamd_application_loop_graph_aware",
        "pygamd_application_loop_is_primary", "pygamd_application_loop_is_eager_debug",
        "cuda_graph_mode", "cuda_graph_strict", "cuda_graph_strict_min_coverage",
        "cuda_graph_strict_max_eager_steps",
        "cuda_graph_diag_enabled", "cuda_graph_diag_env_value",
        "cuda_graph_disabled_internal", "cuda_graph_full_capture_count",
        "cuda_graph_full_replay_count", "cuda_graph_tail_capture_count",
        "cuda_graph_tail_replay_count", "cuda_graph_full_step_count",
        "cuda_graph_tail_step_count", "cuda_graph_eager_step_count",
        "cuda_graph_graph_step_count", "cuda_graph_total_policy_step_count",
        "cuda_graph_coverage_pct", "cuda_graph_fallback_count",
        "cuda_graph_fallback_reasons",
    ):
        if nm in _meta and nm not in snap:
            snap[nm] = _meta.get(nm)
    # coarse run-level timing of the Python step loop (set by application.run)
    snap["python_step_loop_wall_s"] = _times.get("python_step_loop_wall_s", None)
    # conservative kernel-launch estimate = sum of per-step launch-ish counters
    est = sum(_counters.get(k, 0) for k in _LAUNCH_COUNTERS)
    snap["kernel_launch_estimate"] = est if any(k in _counters for k in _LAUNCH_COUNTERS) else None
    # neighbor-list diagnostics: surface meta diagnostic fields at top level (null if not set).
    for nm in _META_FIELDS:
        snap[nm] = _meta.get(nm, None)
    # neighbor-list diagnostics: rebuild ratio = rebuild / check.  >0.8 => list rebuilds almost
    # every check (the DPD signature).  Null when no checks ran.
    chk = _counters.get("nlist_check_count", 0)
    reb = _counters.get("nlist_rebuild_count", 0)
    snap["nlist_rebuild_ratio"] = (reb / chk) if chk else None
    # engine timing: merge the opt-in timing decomposition.  wall_timed_region_ms is
    # the run-level wall of app.run's measured region (register + begin force/sort +
    # dumps + main step loop; set by application.run) so the timing layer can
    # compute unattributed_ms = wall - sum(profiled) over the SAME window.
    if _timing.ENABLED:
        loop_wall = _times.get("python_step_loop_wall_s", None)
        wall_ms = (loop_wall * 1000.0) if isinstance(loop_wall, (int, float)) else None
        snap["timing"] = _timing.snapshot(wall_timed_region_ms=wall_ms)
    # launch trace: merge the opt-in launch / step-loop trace (counts only, no ms).
    if _ltrace.LAUNCH:
        snap["launch_trace"] = _ltrace.snapshot()
    return snap


def dump_json(path, extra=None):
    '''Write the engine metrics profile JSON (counters + run metadata) to `path`.'''
    rec = snapshot()
    rec["profile_enabled"] = ENABLED
    rec["wall_clock_iso"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    if extra:
        rec.update(extra)
    with open(path, "w") as f:
        json.dump(rec, f, indent=2)
    return rec
