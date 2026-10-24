'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License.
	See the General Public License v3 for more details.
CORRESPONDENCE
	Dr. You-Liang Zhu  Email: ylzhu@pygamd.com
'''

# -----------------------------------------------------------------------------
# CuPy-native neighbor list.  Pure cupy.RawKernel; NO numba.  Main line:
#   CuInfo -> CuPy clist -> CuPy nlist -> CuPy RawPair -> CuPy RawNVT
# The legacy Numba implementation moved to archive_legacy/pygamd/plists/nlist.py.
#
# Builds a FULL neighbor list (each i lists every j within rcut+rbuff) using the
# CuPy cell list, with orthorhombic PBC minimum image.  Optional bond / angle /
# dihedral exclusions are tag-space filters applied during neighbor build.  A
# periodic rbuff/2 displacement check decides when to rebuild; neighbor capacity
# overflow grows the per-particle width and rebuilds.  Outputs are CuPy int32
# arrays compatible with RawPair: d_neighbor_list (npa, width) and
# d_neighbor_size (npa,).
#
# Control surface kept for application / CUDA Graph / RawPair:
#   calculate(build_list=False, timestep=None)
#   will_do_gpu_check(build_list=False)  -> read-only prediction
#   did_gpu_check_this_call, did_rebuild_this_call
#   rebuild_count, check_count, last_build_timestep
#   d_neighbor_list, d_neighbor_size
#
# DEFERRED / not yet supported:
#   Future work: further-check, rebuild_on_check experiment, triclinic box, multi-GPU.
#   The check period is configurable via PYGAMD_NLIST_CHECK_PERIOD (default 1).
#
# cupy is imported lazily (kernels compiled on first build); importing this
# module never requires a CUDA driver.
# -----------------------------------------------------------------------------

import os
import math
import numpy as np

import pygamd.plists.clist as clist
from pygamd.backend.cupy import import_cupy
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace


_DEFAULT_CHECK_PERIOD = 1


def _read_check_period():
	'''Resolve PYGAMD_NLIST_CHECK_PERIOD robustly.

	Returns (period, note).  Unset -> default (1), no note.  Any invalid value
	(empty string, 0, negative, non-integer) -> default (1) with a diagnostic
	note string (only surfaced under PYGAMD_ENGINE_METRICS).  This NEVER raises
	and NEVER changes the physics/integrator -- it only changes the displacement
	check interval.  GALAMOST / numba are unaffected (they never read this env).'''
	raw = os.environ.get("PYGAMD_NLIST_CHECK_PERIOD", None)
	if raw is None:
		return _DEFAULT_CHECK_PERIOD, None
	try:
		period = int(str(raw).strip())
	except (TypeError, ValueError):
		return (_DEFAULT_CHECK_PERIOD,
		        "PYGAMD_NLIST_CHECK_PERIOD=%r is not an integer; "
		        "falling back to default %d" % (raw, _DEFAULT_CHECK_PERIOD))
	if period < 1:
		return (_DEFAULT_CHECK_PERIOD,
		        "PYGAMD_NLIST_CHECK_PERIOD=%r is < 1; "
		        "falling back to default %d" % (raw, _DEFAULT_CHECK_PERIOD))
	return period, None


_DEFAULT_RBUFF_SCALE = 0.1   # rbuff = rcut * 0.1 (unchanged default)
_MAX_RBUFF_SCALE = 1.0       # opt-in safe upper bound; above this -> fall back


def _read_rbuff_scale():
	'''Resolve PYGAMD_NLIST_RBUFF_SCALE robustly (neighbor-list skin policy, opt-in).

	Returns (scale, source, note):
	  scale  -- float to use for rbuff = rcut * scale, or None to KEEP the
	            caller-provided rbuff for the unset/default path.
	  source -- "default" (unset), "env" (valid), or "fallback" (invalid).
	  note   -- diagnostic string for invalid input (surfaced only under
	            PYGAMD_ENGINE_METRICS), else None.

	Unset -> (None, "default", None) so the default path is byte-identical.
	Valid 0 < scale <= 1.0 -> (scale, "env", None).
	Invalid (empty / non-number / NaN / <= 0 / > 1.0) ->
	(_DEFAULT_RBUFF_SCALE, "fallback", note), explicitly restoring rcut*0.1 even
	if a future caller supplies a different rbuff.  NEVER raises.  Only affects CuPy pygamd; GALAMOST /
	numba never read this variable, and no default default runner sets it.'''
	raw = os.environ.get("PYGAMD_NLIST_RBUFF_SCALE", None)
	if raw is None:
		return None, "default", None
	try:
		scale = float(str(raw).strip())
	except (TypeError, ValueError):
		return (_DEFAULT_RBUFF_SCALE, "fallback",
		        "PYGAMD_NLIST_RBUFF_SCALE=%r is not a number; falling back to "
		        "default %.3f" % (raw, _DEFAULT_RBUFF_SCALE))
	# NaN fails every comparison, so `scale > 0.0` is False -> fallback.
	if not (scale > 0.0):
		return (_DEFAULT_RBUFF_SCALE, "fallback",
		        "PYGAMD_NLIST_RBUFF_SCALE=%r is <= 0 (or NaN); falling back to "
		        "default %.3f" % (raw, _DEFAULT_RBUFF_SCALE))
	if scale > _MAX_RBUFF_SCALE:
		return (_DEFAULT_RBUFF_SCALE, "fallback",
		        "PYGAMD_NLIST_RBUFF_SCALE=%r exceeds safe max %.2f; falling back "
		        "to default %.3f" % (raw, _MAX_RBUFF_SCALE, _DEFAULT_RBUFF_SCALE))
	return scale, "env", None


def _read_fused_status():
	'''Resolve PYGAMD_NLIST_FUSED_STATUS (fused-status, opt-in).  Returns a bool.

	When ON, the per-rebuild clist-situation read and the neighbor-build width
	read are fused into a SINGLE device->host copy (2 syncs -> 1) in the common
	no-overflow case, with a safe regrow+redo fallback that keeps every clist /
	width safety check.  Default (unset / 0 / false / off) -> OFF, so the default
	path is byte-identical.  Only affects CuPy pygamd; GALAMOST / numba never read
	it, and no default default runner sets it.'''
	v = os.environ.get("PYGAMD_NLIST_FUSED_STATUS", "")
	return v not in ("", "0", "false", "False", "no", "off", "NO", "OFF")


# nlist_check : situation[0]=1 if any particle moved > rbuff/2 since last build.
# nlist_build : situation[1]=max(needed width) when the neighbor buffer overflows.
_NLIST_SRC = r'''
extern "C" __global__ void nlist_check(
		const int   npa,
		const float* __restrict__ pos,
		const float* __restrict__ box,
		const float* __restrict__ last_pos,
		const float  rbuff_half_rsq,
		int*         __restrict__ situation)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float dx = pos[i * 4 + 0] - last_pos[i * 3 + 0];
		float dy = pos[i * 4 + 1] - last_pos[i * 3 + 1];
		float dz = pos[i * 4 + 2] - last_pos[i * 3 + 2];
		float bx = box[0], by = box[1], bz = box[2];
		if      (dx >=  bx * 0.5f) dx -= bx; else if (dx < -bx * 0.5f) dx += bx;
		if      (dy >=  by * 0.5f) dy -= by; else if (dy < -by * 0.5f) dy += by;
		if      (dz >=  bz * 0.5f) dz -= bz; else if (dz < -bz * 0.5f) dz += bz;
		float rsq = dx * dx + dy * dy + dz * dz;
		if (rsq > rbuff_half_rsq) situation[0] = 1;
	}
}

extern "C" __global__ void nlist_build(
		const int   npa,
		const float* __restrict__ pos,
		const int*   __restrict__ dim,
		const float* __restrict__ box_low,
		const float* __restrict__ inv_width,
		const float* __restrict__ box,
		const int*   __restrict__ cell_adj,
		const int*   __restrict__ cell_size,
		const int*   __restrict__ cell_members,
		const int    nmax_cell,
		const int*   __restrict__ tag,
		const int*   __restrict__ excl_size,
		const int*   __restrict__ excl_list,
		const int    excl_width,
		int*         __restrict__ neighbor_size,
		int*         __restrict__ neighbor_list,
		float*       __restrict__ last_pos,
		const float  cutoff_rsq,
		const int    nlist_width,
		int*         __restrict__ situation)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float pix = pos[i * 4 + 0];
		float piy = pos[i * 4 + 1];
		float piz = pos[i * 4 + 2];
		int tagi = tag[i];

		float lx = box_low[0], ly = box_low[1], lz = box_low[2];
		// In the opt-in fused path (PYGAMD_NLIST_FUSED_STATUS) this kernel may run
		// BEFORE the host has read the clist status and raised on a NaN / out-of-box
		// particle, so it must defend against that transient bad state: an invalid
		// pix/piy/piz would otherwise produce an out-of-range cell_id and an
		// out-of-bounds read of cell_adj.  These guards only skip the offending
		// particle (no neighbors) for the discarded intermediate list; the REAL
		// NaN / out-of-box error is still raised by clist.resolve_situation_host()
		// from fused_status[0:2] right after.  In the default path clist.calculate()
		// already host-synced and raised before this kernel runs, so positions are
		// always valid here and these guards never trigger (no behavior change).
		if (isnan(pix) || isnan(piy) || isnan(piz))
		{
			neighbor_size[i] = 0;
			return;
		}
		if (pix < lx || pix >= -lx || piy < ly || piy >= -ly || piz < lz || piz >= -lz)
		{
			neighbor_size[i] = 0;
			return;
		}
		int ix = (int)((pix - lx) * inv_width[0]);
		int iy = (int)((piy - ly) * inv_width[1]);
		int iz = (int)((piz - lz) * inv_width[2]);
		if (ix == dim[0]) ix = 0;
		if (iy == dim[1]) iy = 0;
		if (iz == dim[2]) iz = 0;
		// Final clamp guard: even after the above checks, float rounding could in
		// principle leave an index out of [0, dim); avoid an OOB cell_adj read.
		if (ix < 0 || ix >= dim[0] || iy < 0 || iy >= dim[1] || iz < 0 || iz >= dim[2])
		{
			neighbor_size[i] = 0;
			return;
		}
		int cell_id = iz + dim[2] * (iy + ix * dim[1]);

		float bx = box[0], by = box[1], bz = box[2];
		float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;

		int nneigh = 0;
		int needed = 0;
		for (int ci = 0; ci < 27; ci++)
		{
			int ncell_id = cell_adj[cell_id * 27 + ci];
			if (ncell_id < 0) continue;
			int csize = cell_size[ncell_id];
			// Clamp to the per-cell capacity so we never read past cell_members.
			// In the default path clist overflow is already resolved before this
			// kernel runs, so cell_size <= nmax_cell and this clamp is a no-op.
			// In the opt-in fused path (PYGAMD_NLIST_FUSED_STATUS) this kernel may
			// run BEFORE the host sees a clist overflow, so cell_size can exceed
			// nmax_cell; clamping prevents an out-of-bounds read.  Correctness is
			// preserved because the host still detects the overflow via
			// fused_status[2] (max cell occupancy), grows clist.d_cell_members, and
			// redoes the build -- this truncated intermediate nlist is discarded
			// and never reaches the force kernels.
			if (csize > nmax_cell) csize = nmax_cell;
			for (int s = 0; s < csize; s++)
			{
				int j = cell_members[ncell_id * nmax_cell + s];
				if (j == i) continue;
				if (excl_width > 0)
				{
					int tagj = tag[j];
					int es = excl_size[tagi];
					bool excluded = false;
					for (int ei = 0; ei < es; ei++)
					{
						if (excl_list[tagi * excl_width + ei] == tagj)
						{
							excluded = true;
							break;
						}
					}
					if (excluded) continue;
				}
				float dx = pix - pos[j * 4 + 0];
				float dy = piy - pos[j * 4 + 1];
				float dz = piz - pos[j * 4 + 2];
				if      (dx >=  bxh) dx -= bx; else if (dx < -bxh) dx += bx;
				if      (dy >=  byh) dy -= by; else if (dy < -byh) dy += by;
				if      (dz >=  bzh) dz -= bz; else if (dz < -bzh) dz += bz;
				float rsq = dx * dx + dy * dy + dz * dz;
				if (rsq < cutoff_rsq)
				{
					if (nneigh < nlist_width) neighbor_list[i * nlist_width + nneigh] = j;
					else                      needed = nneigh + 1;
					nneigh++;
				}
			}
		}
		neighbor_size[i] = nneigh;
		last_pos[i * 3 + 0] = pix;
		last_pos[i * 3 + 1] = piy;
		last_pos[i * 3 + 2] = piz;
		if (needed > 0) atomicMax(&situation[1], needed);
	}
}
'''


# OPT-IN profiling variant of nlist_check: identical displacement test, but ALSO
# records the maximum displacement-squared (scaled to int) in situation[2] via a
# single atomicMax.  Compiled and used ONLY when PYGAMD_ENGINE_METRICS is on, so
# the default (non-profile) path keeps the lean kernel above byte-for-byte.  The
# value is read from the SAME host copy of `situation` we already fetch for the
# displacement flag, so it adds NO extra host sync.  Scale = 1e6 (rsq is small:
# threshold rsq = (rbuff/2)^2); max displacement = sqrt(situation[2] / 1e6).
_NLIST_CHECK_PROF_SRC = r'''
extern "C" __global__ void nlist_check_prof(
		const int   npa,
		const float* __restrict__ pos,
		const float* __restrict__ box,
		const float* __restrict__ last_pos,
		const float  rbuff_half_rsq,
		int*         __restrict__ situation)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float dx = pos[i * 4 + 0] - last_pos[i * 3 + 0];
		float dy = pos[i * 4 + 1] - last_pos[i * 3 + 1];
		float dz = pos[i * 4 + 2] - last_pos[i * 3 + 2];
		float bx = box[0], by = box[1], bz = box[2];
		if      (dx >=  bx * 0.5f) dx -= bx; else if (dx < -bx * 0.5f) dx += bx;
		if      (dy >=  by * 0.5f) dy -= by; else if (dy < -by * 0.5f) dy += by;
		if      (dz >=  bz * 0.5f) dz -= bz; else if (dz < -bz * 0.5f) dz += bz;
		float rsq = dx * dx + dy * dy + dz * dz;
		if (rsq > rbuff_half_rsq) situation[0] = 1;
		atomicMax(&situation[2], (int)(rsq * 1.0e6f));
	}
}
'''
_NLIST_DISP_SCALE = 1.0e6   # mirror of the scale in _NLIST_CHECK_PROF_SRC


class nlist:
	def __init__(self, info, rcut, rbuff, exclusion):
		self.info = info
		self.cu = getattr(info, "cu", None)
		self.rcut = rcut
		self.exclusion = tuple(exclusion) if exclusion else None
		# Opt-in rbuff/skin policy (neighbor-list skin policy): rbuff = rcut * PYGAMD_NLIST_RBUFF_SCALE.
		# Unset -> keep the caller-provided rbuff (rcut*0.1) exactly, so the default
		# path is byte-identical.  Invalid env -> explicitly fall back to rcut*0.1.
		_scale, self._rbuff_scale_source, _scale_note = _read_rbuff_scale()
		if _scale is not None and rcut > 0:
			rbuff = float(rcut) * _scale
		self._rbuff_scale_note = _scale_note
		self.rbuff = rbuff
		# effective scale derived from the final rbuff (covers default/env/fallback)
		self._rbuff_scale_effective = (float(rbuff) / float(rcut)) if rcut > 0 else None
		self.cutoff_rsq = float((rcut + rbuff) * (rcut + rbuff))
		self.rbuff_half_rsq = float((rbuff / 2.0) * (rbuff / 2.0))
		self.block_size = 64
		self.nblocks = math.ceil(self.info.npa / self.block_size)

		# CuPy cell list sized to the search cutoff (rcut + rbuff)
		self.clist = clist.clist(info, rcut + rbuff)

		self.cp = None
		self._k_check = None
		self._k_check_prof = None           # opt-in profiling check kernel (ENABLED only)
		self._k_build = None
		self.d_excl_size = None
		self.d_excl_list = None
		self.excl_width = 0
		self.nmax = 8                       # per-particle neighbor capacity
		self.d_neighbor_list = None         # (npa, nmax) int32
		self.d_neighbor_size = None         # (npa,) int32
		self.d_last_pos = None              # (npa, 3) float32
		self.d_situation = None             # (3,) int32
		# fused-status opt-in fused build-status: clist[0:3] + nbuild[3:6] in ONE buffer.
		self.fused_status = _read_fused_status()
		self.d_fused_status = None          # (6,) int32, allocated only when fused

		# periodic-check state + read-only flags / counters.  Invalid env values
		# fall back to the default and are noted (only under PYGAMD_ENGINE_METRICS).
		self.check_period, _cp_note = _read_check_period()
		if _sp.ENABLED:
			_sp.set_meta("nlist_check_period_effective", int(self.check_period))
			_sp.set_meta("nlist_skin_threshold", float(self.rbuff / 2.0))
			_sp.set_meta("nlist_rbuff_scale_effective", self._rbuff_scale_effective)
			_sp.set_meta("nlist_rbuff_scale_source", self._rbuff_scale_source)
			_sp.set_meta("nlist_fused_status_enabled", bool(self.fused_status))
			if _cp_note:
				_sp.note(_cp_note)
			if self._rbuff_scale_note:
				_sp.note(self._rbuff_scale_note)
			if self._rbuff_scale_source == "env":
				_sp.note("PYGAMD_NLIST_RBUFF_SCALE active (rbuff=%.4g, scale=%.4g): "
				         "opt-in larger skin; verify runtime/metrics before adopting"
				         % (float(self.rbuff), float(self._rbuff_scale_effective)))
			if self.fused_status:
				_sp.note("PYGAMD_NLIST_FUSED_STATUS active: clist+width status reads fused into "
				         "one host copy per rebuild (2 syncs -> 1); safe regrow+redo fallback kept")
		self._built = False
		self._last_check_call = -1
		self._last_box = None
		self.call = 0
		self.did_gpu_check_this_call = False
		self.did_rebuild_this_call = False
		self.rebuild_count = 0
		self.check_count = 0
		self.last_build_timestep = -1

	# --- exclusions -----------------------------------------------------------
	@staticmethod
	def _add_excl_pair(excl_sets, a, b):
		if int(a) == int(b):
			return
		excl_sets[int(a)].add(int(b))
		excl_sets[int(b)].add(int(a))

	def _collect_exclusion_pairs_host(self):
		'''Build a symmetric tag-space exclusion table on host.

		The table is generated once from CuInfo topology tags and then consumed on
		GPU as tag -> excluded neighbor tags.  Because it is tag-space, Sort does
		not require any exclusion-table refresh; the nlist kernel reads d_tag[i]
		and d_tag[j] and filters by tag identity.
		'''
		npa = int(self.info.npa)
		excl_sets = [set() for _ in range(npa)]
		kinds = set(self.exclusion or ())
		cu = getattr(self.info, "cu", None)
		topology = getattr(cu, "topology", {}) if cu is not None else {}

		if "bond" in kinds and "bond" in topology:
			for a, b in topology["bond"].host_tags[:, :2]:
				self._add_excl_pair(excl_sets, a, b)
		if "angle" in kinds and "angle" in topology:
			# Angle exclusion means the 1-3 endpoint pair.  Bonded 1-2 pairs are
			# only excluded when "bond" is also requested.
			for a, _b, c in topology["angle"].host_tags[:, :3]:
				self._add_excl_pair(excl_sets, a, c)
		if "dihedral" in kinds and "dihedral" in topology:
			# Dihedral exclusion means the 1-4 endpoint pair.
			for a, _b, _c, d in topology["dihedral"].host_tags[:, :4]:
				self._add_excl_pair(excl_sets, a, d)

		max_width = max((len(x) for x in excl_sets), default=0)
		sizes = np.asarray([len(x) for x in excl_sets], dtype=np.int32)
		if max_width > 0:
			lst = np.zeros((npa, max_width), dtype=np.int32)
			for tag_id, vals in enumerate(excl_sets):
				if vals:
					lst[tag_id, :len(vals)] = np.asarray(sorted(vals), dtype=np.int32)
		else:
			lst = np.zeros((npa, 1), dtype=np.int32)   # dummy array for kernel args
		return sizes, lst, int(max_width)


	# --- device setup ---------------------------------------------------------
	def _ensure_init(self):
		if self.cp is not None:
			return
		cp = import_cupy("CuPy neighbor list")
		self.cp = cp
		mod = cp.RawModule(code=_NLIST_SRC)
		self._k_check = mod.get_function("nlist_check")
		self._k_build = mod.get_function("nlist_build")
		# Opt-in: compile the displacement-recording check kernel ONLY when engine metrics
		# profiling is enabled, so the default benchmark path is unchanged.
		if _sp.ENABLED:
			prof_mod = cp.RawModule(code=_NLIST_CHECK_PROF_SRC)
			self._k_check_prof = prof_mod.get_function("nlist_check_prof")
		npa = int(self.info.npa)
		self.d_neighbor_list = cp.zeros((npa, self.nmax), dtype=cp.int32)
		self.d_neighbor_size = cp.zeros(npa, dtype=cp.int32)
		self.d_last_pos = cp.zeros((npa, 3), dtype=cp.float32)
		self.d_situation = cp.zeros(3, dtype=cp.int32)
		excl_size, excl_list, excl_width = self._collect_exclusion_pairs_host()
		self.d_excl_size = cp.asarray(excl_size, dtype=cp.int32)
		self.d_excl_list = cp.asarray(excl_list, dtype=cp.int32)
		self.excl_width = int(excl_width)
		if _sp.ENABLED:
			_sp.set_meta("nlist_exclusion_kinds", ",".join(self.exclusion or ()))
			_sp.set_meta("nlist_exclusion_max_per_tag", int(self.excl_width))
			_sp.set_meta("nlist_exclusion_total_directed", int(excl_size.sum()))
		# fused-status: fused status buffer [clist(0:3) | neighbor-build(3:6)].
		if self.fused_status:
			self.d_fused_status = cp.zeros(6, dtype=cp.int32)

	def _dev(self, name, fallback_attr):
		if self.cu is not None:
			return self.cu.dev(name)
		from pygamd.backend.cupy import as_cupy_array
		return as_cupy_array(getattr(self.info, fallback_attr), ravel=True)

	# --- read-only check prediction (CUDA graph) ------------------------------
	def will_do_gpu_check(self, build_list=False):
		'''Predict whether the NEXT calculate(build_list) will run a GPU check
		(i.e. NOT a pure periodic skip).  Mutates nothing.'''
		if build_list:
			return True
		if not self._built:
			return True
		b = self.info.box
		if (self._last_box is None or b[0] != self._last_box[0]
		        or b[1] != self._last_box[1] or b[2] != self._last_box[2]):
			return True
		call_next = self.call + 1
		return (call_next - self._last_check_call) >= self.check_period

	def _decide_check(self):
		if not self._built:
			return True
		b = self.info.box
		if (self._last_box is None or b[0] != self._last_box[0]
		        or b[1] != self._last_box[1] or b[2] != self._last_box[2]):
			return True
		return (self.call - self._last_check_call) >= self.check_period

	# --- main entry -----------------------------------------------------------
	def calculate(self, build_list=False, timestep=None):
		prof.begin("nlist_calculate_total")
		prof.count("nlist_calculate_total")
		self._ensure_init()
		self.call += 1
		self.did_gpu_check_this_call = False
		self.did_rebuild_this_call = False

		if not (build_list or self._decide_check()):
			prof.count("nlist_periodic_skip")
			prof.end("nlist_calculate_total")
			return

		# committing to a GPU check this call
		self.did_gpu_check_this_call = True
		self.check_count += 1
		if _sp.ENABLED:
			_sp.inc("nlist_check_count")
			_sp.set_meta("nlist_last_check_step",
			             int(timestep) if timestep is not None else int(self.call))
		prof.count("nlist_periodic_check")
		self._last_check_call = self.call
		b = self.info.box
		self._last_box = (b[0], b[1], b[2])

		# rebuild reason classification (diagnostic only; does not change behavior):
		#   forced            -> caller requested build_list=True
		#   periodic_check    -> first build / structural (not yet built)
		#   displacement      -> a particle moved > rbuff/2 since last build
		#   width_realloc     -> extra build iteration due to capacity overflow
		need_rebuild = bool(build_list) or (not self._built)
		rebuild_reason = ("forced" if build_list else
		                  ("periodic_check" if not self._built else None))
		if not need_rebuild:
			cp = self.cp
			self.d_situation.fill(0)
			use_prof = bool(_sp.ENABLED and self._k_check_prof is not None)
			kern = self._k_check_prof if use_prof else self._k_check
			prof.begin("nlist_check_build_kernel")
			if _lt.LAUNCH: _lt.launch("nlist_check_kernel")
			if _t.ENABLED: _t.gpu_begin("nlist_check_kernel")
			kern((self.nblocks,), (self.block_size,),
			     (np.int32(self.info.npa), self._dev("pos", "d_pos"),
			      self._dev("box", "d_box"), self.d_last_pos.ravel(),
			      np.float32(self.rbuff_half_rsq), self.d_situation))
			if _t.ENABLED: _t.gpu_end("nlist_check_kernel")
			prof.end("nlist_check_build_kernel")
			prof.begin("nlist_check_sync")
			if _t.ENABLED: _t.host_begin("nlist_displacement_copy")
			sit = cp.asnumpy(self.d_situation)        # single device->host copy
			if _t.ENABLED: _t.host_end("nlist_displacement_copy")
			flag = int(sit[0])
			if _sp.ENABLED:
				_sp.inc("host_sync_count")            # device->host check sync
				_sp.inc("host_sync_reason_displacement_check")
				if use_prof:
					# max displacement is read from the SAME host copy: no extra sync.
					maxrsq = float(sit[2]) / _NLIST_DISP_SCALE
					if maxrsq > 0.0:
						_sp.set_meta_max("nlist_max_displacement_seen", math.sqrt(maxrsq))
			prof.end("nlist_check_sync")
			if flag == 1:
				need_rebuild = True
				rebuild_reason = "displacement"

		if need_rebuild:
			self._rebuild(timestep, rebuild_reason)

		self._built = True
		prof.end("nlist_calculate_total")

	def _rebuild(self, timestep, reason="forced"):
		cp = self.cp
		self.did_rebuild_this_call = True
		self.rebuild_count += 1
		if _sp.ENABLED:
			_sp.inc("nlist_rebuild_count")
			_sp.inc("nlist_rebuild_due_to_%s" % (reason or "periodic_check"))
			_sp.set_meta("nlist_last_rebuild_step",
			             int(timestep) if timestep is not None else int(self.call))
		self.last_build_timestep = int(timestep) if timestep is not None else int(self.call)
		prof.count("nlist_rebuild_count")

		building = True
		while building:
			if self.fused_status and self.d_fused_status is not None:
				building = self._build_once_fused(timestep)
				continue
			self.clist.calculate()
			self.d_situation.fill(0)
			nlist_width = int(self.d_neighbor_list.shape[1])
			prof.begin("nlist_neighbor_build_kernel")
			if _lt.LAUNCH: _lt.launch("nlist_build_kernel")
			if _t.ENABLED: _t.gpu_begin("nlist_build_kernel")
			self._k_build((self.nblocks,), (self.block_size,),
			              (np.int32(self.info.npa), self._dev("pos", "d_pos"),
			               self.clist.d_dim, self.clist.d_box_low_boundary,
			               self.clist.d_inv_width, self._dev("box", "d_box"),
			               self.clist.d_cell_adj, self.clist.d_cell_size,
			               self.clist.d_cell_members, np.int32(self.clist.nmax),
			               self._dev("tag", "d_tag"), self.d_excl_size, self.d_excl_list,
			               np.int32(self.excl_width), self.d_neighbor_size, self.d_neighbor_list,
			               self.d_last_pos.ravel(), np.float32(self.cutoff_rsq),
			               np.int32(nlist_width), self.d_situation))
			if _t.ENABLED: _t.gpu_end("nlist_build_kernel")
			prof.end("nlist_neighbor_build_kernel")
			prof.begin("nlist_neighbor_build_sync")
			if _t.ENABLED: _t.host_begin("nlist_width_copy")
			needed = int(cp.asnumpy(self.d_situation)[1])
			if _t.ENABLED: _t.host_end("nlist_width_copy")
			if _sp.ENABLED:
				_sp.inc("host_sync_count")            # device->host neighbor-build sync
				_sp.inc("host_sync_reason_width_check")
			prof.end("nlist_neighbor_build_sync")
			if needed > self.nmax:
				prof.count("nlist_expand_count")
				if _sp.ENABLED:
					_sp.inc("nlist_width_realloc_count")
					# an overflow forces an EXTRA build iteration -> width_realloc reason
					_sp.inc("nlist_rebuild_due_to_width_realloc")
				self.nmax = int(needed + 8 - (needed & 7))
				if _t.ENABLED: _t.host_begin("nlist_realloc")
				self.d_neighbor_list = cp.zeros((int(self.info.npa), self.nmax), dtype=cp.int32)
				if _t.ENABLED: _t.host_end("nlist_realloc")
			else:
				building = False

	def _build_once_fused(self, timestep):
		'''fused-status opt-in fused build-status iteration.

		Runs the clist cell-build kernel (status into d_fused_status[0:3]) and the
		neighbor-build kernel (needed-width into d_fused_status[3:6]) back-to-back
		on the same stream, then reads BOTH with a SINGLE device->host copy
		(2 syncs -> 1).  Correctness is preserved: clist NaN / out-of-box / cell
		overflow and per-particle width overflow are still checked from that copy,
		and on any overflow we regrow + return True so the caller redoes the build
		(the truncated intermediate list is discarded before forces).  Returns True
		if another iteration is needed (overflow), else False.'''
		cp = self.cp
		self.d_fused_status.fill(0)
		# clist writes [nan, oob, maxocc] into the first 3 slots (no internal sync).
		# Because we do NOT sync here, a clist cell overflow is not yet visible on
		# the host; nlist_build runs optimistically.  The nlist_build kernel clamps
		# each cell's csize to nmax_cell so an overflowed cell cannot cause an
		# out-of-bounds read of cell_members.  The overflow is then detected from
		# the fused copy (fused[2]) below, which regrows clist and forces a redo.
		self.clist.build_into(self.d_fused_status[0:3])
		nlist_width = int(self.d_neighbor_list.shape[1])
		prof.begin("nlist_neighbor_build_kernel")
		if _lt.LAUNCH: _lt.launch("nlist_build_kernel")
		if _t.ENABLED: _t.gpu_begin("nlist_build_kernel")
		self._k_build((self.nblocks,), (self.block_size,),
		              (np.int32(self.info.npa), self._dev("pos", "d_pos"),
		               self.clist.d_dim, self.clist.d_box_low_boundary,
		               self.clist.d_inv_width, self._dev("box", "d_box"),
		               self.clist.d_cell_adj, self.clist.d_cell_size,
		               self.clist.d_cell_members, np.int32(self.clist.nmax),
		               self._dev("tag", "d_tag"), self.d_excl_size, self.d_excl_list,
		               np.int32(self.excl_width), self.d_neighbor_size, self.d_neighbor_list,
		               self.d_last_pos.ravel(), np.float32(self.cutoff_rsq),
		               np.int32(nlist_width), self.d_fused_status[3:6]))
		if _t.ENABLED: _t.gpu_end("nlist_build_kernel")
		prof.end("nlist_neighbor_build_kernel")
		prof.begin("nlist_fused_status_sync")
		if _t.ENABLED: _t.host_begin("nlist_fused_status_copy")
		fused = cp.asnumpy(self.d_fused_status)        # ONE copy for clist + width
		if _t.ENABLED: _t.host_end("nlist_fused_status_copy")
		if _sp.ENABLED:
			_sp.inc("host_sync_count")
			_sp.inc("host_sync_reason_fused_status")
			_sp.inc("fused_status_copy_count")
		prof.end("nlist_fused_status_sync")
		# clist safety: raises on NaN/out-of-box; grows + returns True on overflow
		clist_grew = self.clist.resolve_situation_host(fused[0:3])
		needed = int(fused[4])                          # nbuild situation[1] -> fused[4]
		width_grew = False
		if needed > self.nmax:
			prof.count("nlist_expand_count")
			if _sp.ENABLED:
				_sp.inc("nlist_width_realloc_count")
				_sp.inc("nlist_rebuild_due_to_width_realloc")
			self.nmax = int(needed + 8 - (needed & 7))
			if _t.ENABLED: _t.host_begin("nlist_realloc")
			self.d_neighbor_list = cp.zeros((int(self.info.npa), self.nmax), dtype=cp.int32)
			if _t.ENABLED: _t.host_end("nlist_realloc")
			width_grew = True
		return bool(clist_grew or width_grew)

	def speak(self):
		print("CuPy nlist: npa=%d width=%d rebuilds=%d checks=%d last_build_ts=%d exclusions=%s excl_width=%d"
		      % (int(self.info.npa), self.nmax, self.rebuild_count,
		         self.check_count, self.last_build_timestep,
		         ",".join(self.exclusion or ()), int(self.excl_width)))
