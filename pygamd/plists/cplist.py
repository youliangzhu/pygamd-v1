'''
CuPy-native cluster pair list.

For large systems this builds the same public neighbor-list interface as nlist
(`d_neighbor_list`, `d_neighbor_size`) but traverses candidate partners as
fixed-size particle clusters.  The force backends can therefore keep their
existing kernels while the expensive list rebuild path gets a more GPU-friendly
cluster-pair traversal.

Cluster traversal:
  * build the normal CuPy cell list;
  * split every cell's particle slots into fixed clusters of 4/8/16 particles;
  * launch one CUDA block per source cluster;
  * let the block's cluster-size^2 lanes scan target clusters from the 27
    neighboring cells, writing directed FULL neighbor entries with atomicAdd.

This keeps the nlist control surface intact for CUDA Graph and force backends:
  calculate(build_list=False, timestep=None)
  will_do_gpu_check(build_list=False)
  did_gpu_check_this_call, did_rebuild_this_call
  rebuild_count, check_count, last_build_timestep
  d_neighbor_list, d_neighbor_size
'''

import math
import os
import numpy as np

import pygamd.plists.clist as clist
from pygamd.backend.cupy import import_cupy
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp
from pygamd import engine_timing as _t
from pygamd import launch_trace as _lt
from pygamd.plists.nlist import _read_check_period, _read_rbuff_scale


_MAX_CLUSTER_SIZE = 16
_DEFAULT_CLUSTER_SIZE = 8
_DEFAULT_INITIAL_WIDTH = 64


def _read_cluster_size():
	raw = os.environ.get("PYGAMD_CLUSTER_SIZE", "")
	if raw == "":
		return _DEFAULT_CLUSTER_SIZE
	try:
		size = int(str(raw).strip())
	except (TypeError, ValueError):
		return _DEFAULT_CLUSTER_SIZE
	if size in (4, 8, 16):
		return size
	return _DEFAULT_CLUSTER_SIZE


def _read_initial_width():
	raw = os.environ.get("PYGAMD_CLUSTER_NLIST_INITIAL_WIDTH", "")
	if raw == "":
		return _DEFAULT_INITIAL_WIDTH
	try:
		width = int(str(raw).strip())
	except (TypeError, ValueError):
		return _DEFAULT_INITIAL_WIDTH
	if width < 8:
		return 8
	if width > 4096:
		return 4096
	return int((width + 7) & ~7)


_CPLIST_SRC = r'''
extern "C" __global__ void cplist_check(
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

extern "C" __global__ void cplist_copy_last_pos(
		const int npa,
		const float* __restrict__ pos,
		float* __restrict__ last_pos)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		last_pos[i * 3 + 0] = pos[i * 4 + 0];
		last_pos[i * 3 + 1] = pos[i * 4 + 1];
		last_pos[i * 3 + 2] = pos[i * 4 + 2];
	}
}

extern "C" __global__ void cplist_build(
		const int   ncluster_slots,
		const int   cluster_size,
		const int   clusters_per_cell,
		const float* __restrict__ pos,
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
		const float  cutoff_rsq,
		const int    nlist_width,
		int*         __restrict__ situation)
{
	int source_cluster = blockIdx.x;
	if (source_cluster >= ncluster_slots) return;

	int source_cell = source_cluster / clusters_per_cell;
	int source_slot = source_cluster - source_cell * clusters_per_cell;
	int source_start = source_slot * cluster_size;
	int source_cell_size = cell_size[source_cell];
	int source_count = source_cell_size - source_start;
	if (source_count <= 0) return;
	if (source_count > cluster_size) source_count = cluster_size;

	int tid = threadIdx.x;
	int si_local = tid / cluster_size;
	int tj_local = tid - si_local * cluster_size;
	if (si_local >= source_count) return;

	int i = cell_members[source_cell * nmax_cell + source_start + si_local];
	int tagi = tag[i];
	float pix = pos[i * 4 + 0];
	float piy = pos[i * 4 + 1];
	float piz = pos[i * 4 + 2];

	float bx = box[0], by = box[1], bz = box[2];
	float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;

	for (int ai = 0; ai < 27; ai++)
	{
		int target_cell = cell_adj[source_cell * 27 + ai];
		if (target_cell < 0) continue;
		int target_cell_size = cell_size[target_cell];
		for (int target_slot = 0; target_slot < clusters_per_cell; target_slot++)
		{
			int target_start = target_slot * cluster_size;
			int target_count = target_cell_size - target_start;
			if (target_count <= 0) continue;
			if (target_count > cluster_size) target_count = cluster_size;
			if (tj_local < target_count)
			{
				int j = cell_members[target_cell * nmax_cell + target_start + tj_local];
				if (i != j)
				{
					bool excluded = false;
					if (excl_width > 0)
					{
						int tagj = tag[j];
						int es = excl_size[tagi];
						for (int ei = 0; ei < es; ei++)
						{
							if (excl_list[tagi * excl_width + ei] == tagj)
							{
								excluded = true;
								break;
							}
						}
					}
					if (!excluded)
					{
						float dx = pix - pos[j * 4 + 0];
						float dy = piy - pos[j * 4 + 1];
						float dz = piz - pos[j * 4 + 2];
						if      (dx >=  bxh) dx -= bx; else if (dx < -bxh) dx += bx;
						if      (dy >=  byh) dy -= by; else if (dy < -byh) dy += by;
						if      (dz >=  bzh) dz -= bz; else if (dz < -bzh) dz += bz;
						float rsq = dx * dx + dy * dy + dz * dz;
						if (rsq < cutoff_rsq)
						{
							int slot = atomicAdd(&neighbor_size[i], 1);
							if (slot < nlist_width)
								neighbor_list[i * nlist_width + slot] = j;
							else
								atomicMax(&situation[1], slot + 1);
						}
					}
				}
			}
		}
	}
}
'''


class cplist:
	def __init__(self, info, rcut, rbuff, exclusion):
		self.info = info
		self.cu = getattr(info, "cu", None)
		self.rcut = rcut
		self.exclusion = tuple(exclusion) if exclusion else None
		_scale, self._rbuff_scale_source, _scale_note = _read_rbuff_scale()
		if _scale is not None and rcut > 0:
			rbuff = float(rcut) * _scale
		self._rbuff_scale_note = _scale_note
		self.rbuff = rbuff
		self._rbuff_scale_effective = (float(rbuff) / float(rcut)) if rcut > 0 else None
		self.cutoff_rsq = float((rcut + rbuff) * (rcut + rbuff))
		self.rbuff_half_rsq = float((rbuff / 2.0) * (rbuff / 2.0))

		self.cluster_size = _read_cluster_size()
		if self.cluster_size > _MAX_CLUSTER_SIZE:
			self.cluster_size = _MAX_CLUSTER_SIZE
		self.block_size = int(self.cluster_size * self.cluster_size)
		self.particle_block_size = 128
		self.particle_nblocks = math.ceil(self.info.npa / self.particle_block_size)

		self.clist = clist.clist(info, rcut + rbuff)

		self.cp = None
		self._k_check = None
		self._k_build = None
		self._k_copy_last_pos = None
		self.d_excl_size = None
		self.d_excl_list = None
		self.excl_width = 0
		self.nmax = _read_initial_width()
		self.d_neighbor_list = None
		self.d_neighbor_size = None
		self.d_last_pos = None
		self.d_situation = None
		self.clusters_per_cell = 0
		self.ncluster_slots = 0

		self.check_period, _cp_note = _read_check_period()
		if _sp.ENABLED:
			_sp.set_meta("plist_algorithm", "cluster_pair")
			_sp.set_meta("cluster_pair_cluster_size", int(self.cluster_size))
			_sp.set_meta("cluster_pair_initial_width", int(self.nmax))
			_sp.set_meta("cluster_pair_threshold_active", True)
			_sp.set_meta("nlist_check_period_effective", int(self.check_period))
			_sp.set_meta("nlist_skin_threshold", float(self.rbuff / 2.0))
			_sp.set_meta("nlist_rbuff_scale_effective", self._rbuff_scale_effective)
			_sp.set_meta("nlist_rbuff_scale_source", self._rbuff_scale_source)
			if _cp_note:
				_sp.note(_cp_note)
			if self._rbuff_scale_note:
				_sp.note(self._rbuff_scale_note)

		self._built = False
		self._last_check_call = -1
		self._last_box = None
		self.call = 0
		self.did_gpu_check_this_call = False
		self.did_rebuild_this_call = False
		self.rebuild_count = 0
		self.check_count = 0
		self.last_build_timestep = -1

	@staticmethod
	def _add_excl_pair(excl_sets, a, b):
		if int(a) == int(b):
			return
		excl_sets[int(a)].add(int(b))
		excl_sets[int(b)].add(int(a))

	def _collect_exclusion_pairs_host(self):
		npa = int(self.info.npa)
		excl_sets = [set() for _ in range(npa)]
		kinds = set(self.exclusion or ())
		cu = getattr(self.info, "cu", None)
		topology = getattr(cu, "topology", {}) if cu is not None else {}

		if "bond" in kinds and "bond" in topology:
			for a, b in topology["bond"].host_tags[:, :2]:
				self._add_excl_pair(excl_sets, a, b)
		if "angle" in kinds and "angle" in topology:
			for a, _b, c in topology["angle"].host_tags[:, :3]:
				self._add_excl_pair(excl_sets, a, c)
		if "dihedral" in kinds and "dihedral" in topology:
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
			lst = np.zeros((npa, 1), dtype=np.int32)
		return sizes, lst, int(max_width)

	def _ensure_init(self):
		if self.cp is not None:
			return
		cp = import_cupy("CuPy cluster pair list")
		self.cp = cp
		mod = cp.RawModule(code=_CPLIST_SRC)
		self._k_check = mod.get_function("cplist_check")
		self._k_build = mod.get_function("cplist_build")
		self._k_copy_last_pos = mod.get_function("cplist_copy_last_pos")
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

	def _dev(self, name, fallback_attr):
		if self.cu is not None:
			return self.cu.dev(name)
		from pygamd.backend.cupy import as_cupy_array
		return as_cupy_array(getattr(self.info, fallback_attr), ravel=True)

	def _update_cluster_geometry(self):
		self.clusters_per_cell = int((int(self.clist.nmax) + self.cluster_size - 1) // self.cluster_size)
		if self.clusters_per_cell < 1:
			self.clusters_per_cell = 1
		self.ncluster_slots = int(self.clist.ncell) * int(self.clusters_per_cell)
		if _sp.ENABLED:
			_sp.set_meta("cluster_pair_clusters_per_cell", int(self.clusters_per_cell))
			_sp.set_meta("cluster_pair_cluster_slots", int(self.ncluster_slots))

	def will_do_gpu_check(self, build_list=False):
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

	def calculate(self, build_list=False, timestep=None):
		prof.begin("cplist_calculate_total")
		prof.count("cplist_calculate_total")
		self._ensure_init()
		self.call += 1
		self.did_gpu_check_this_call = False
		self.did_rebuild_this_call = False

		if not (build_list or self._decide_check()):
			prof.count("cplist_periodic_skip")
			prof.end("cplist_calculate_total")
			return

		self.did_gpu_check_this_call = True
		self.check_count += 1
		if _sp.ENABLED:
			_sp.inc("nlist_check_count")
			_sp.inc("cluster_pair_check_count")
			_sp.set_meta("nlist_last_check_step",
			             int(timestep) if timestep is not None else int(self.call))
		prof.count("cplist_periodic_check")
		self._last_check_call = self.call
		b = self.info.box
		self._last_box = (b[0], b[1], b[2])

		need_rebuild = bool(build_list) or (not self._built)
		rebuild_reason = ("forced" if build_list else
		                  ("periodic_check" if not self._built else None))
		if not need_rebuild:
			cp = self.cp
			self.d_situation.fill(0)
			prof.begin("cplist_check_kernel")
			if _lt.LAUNCH: _lt.launch("cplist_check_kernel")
			if _t.ENABLED: _t.gpu_begin("cplist_check_kernel")
			self._k_check((self.particle_nblocks,), (self.particle_block_size,),
			              (np.int32(self.info.npa), self._dev("pos", "d_pos"),
			               self._dev("box", "d_box"), self.d_last_pos.ravel(),
			               np.float32(self.rbuff_half_rsq), self.d_situation))
			if _t.ENABLED: _t.gpu_end("cplist_check_kernel")
			prof.end("cplist_check_kernel")
			prof.begin("cplist_check_sync")
			if _t.ENABLED: _t.host_begin("cplist_displacement_copy")
			sit = cp.asnumpy(self.d_situation)
			if _t.ENABLED: _t.host_end("cplist_displacement_copy")
			if _sp.ENABLED:
				_sp.inc("host_sync_count")
				_sp.inc("host_sync_reason_displacement_check")
			prof.end("cplist_check_sync")
			if int(sit[0]) == 1:
				need_rebuild = True
				rebuild_reason = "displacement"

		if need_rebuild:
			self._rebuild(timestep, rebuild_reason)

		self._built = True
		prof.end("cplist_calculate_total")

	def _rebuild(self, timestep, reason="forced"):
		cp = self.cp
		self.did_rebuild_this_call = True
		self.rebuild_count += 1
		if _sp.ENABLED:
			_sp.inc("nlist_rebuild_count")
			_sp.inc("cluster_pair_rebuild_count")
			_sp.inc("nlist_rebuild_due_to_%s" % (reason or "periodic_check"))
			_sp.set_meta("nlist_last_rebuild_step",
			             int(timestep) if timestep is not None else int(self.call))
		self.last_build_timestep = int(timestep) if timestep is not None else int(self.call)
		prof.count("cplist_rebuild_count")

		self.clist.calculate()
		self._update_cluster_geometry()

		building = True
		while building:
			self.d_neighbor_size.fill(0)
			self.d_situation.fill(0)
			nlist_width = int(self.d_neighbor_list.shape[1])
			prof.begin("cplist_neighbor_build_kernel")
			if _lt.LAUNCH: _lt.launch("cplist_build_kernel")
			if _t.ENABLED: _t.gpu_begin("cplist_build_kernel")
			self._k_build((int(self.ncluster_slots),), (int(self.block_size),),
			              (np.int32(self.ncluster_slots), np.int32(self.cluster_size),
			               np.int32(self.clusters_per_cell), self._dev("pos", "d_pos"),
			               self._dev("box", "d_box"), self.clist.d_cell_adj,
			               self.clist.d_cell_size, self.clist.d_cell_members,
			               np.int32(self.clist.nmax), self._dev("tag", "d_tag"),
			               self.d_excl_size, self.d_excl_list, np.int32(self.excl_width),
			               self.d_neighbor_size, self.d_neighbor_list,
			               np.float32(self.cutoff_rsq), np.int32(nlist_width),
			               self.d_situation))
			if _t.ENABLED: _t.gpu_end("cplist_build_kernel")
			prof.end("cplist_neighbor_build_kernel")
			prof.begin("cplist_neighbor_build_sync")
			if _t.ENABLED: _t.host_begin("cplist_width_copy")
			needed = int(cp.asnumpy(self.d_situation)[1])
			if _t.ENABLED: _t.host_end("cplist_width_copy")
			if _sp.ENABLED:
				_sp.inc("host_sync_count")
				_sp.inc("host_sync_reason_width_check")
			prof.end("cplist_neighbor_build_sync")
			if needed > self.nmax:
				prof.count("cplist_expand_count")
				if _sp.ENABLED:
					_sp.inc("nlist_width_realloc_count")
					_sp.inc("cluster_pair_width_realloc_count")
					_sp.inc("nlist_rebuild_due_to_width_realloc")
				self.nmax = int((needed + 7) & ~7)
				if _t.ENABLED: _t.host_begin("cplist_realloc")
				self.d_neighbor_list = cp.zeros((int(self.info.npa), self.nmax), dtype=cp.int32)
				if _t.ENABLED: _t.host_end("cplist_realloc")
			else:
				building = False

		prof.begin("cplist_copy_last_pos_kernel")
		if _lt.LAUNCH: _lt.launch("cplist_copy_last_pos_kernel")
		if _t.ENABLED: _t.gpu_begin("cplist_copy_last_pos_kernel")
		self._k_copy_last_pos((self.particle_nblocks,), (self.particle_block_size,),
		                      (np.int32(self.info.npa), self._dev("pos", "d_pos"),
		                       self.d_last_pos.ravel()))
		if _t.ENABLED: _t.gpu_end("cplist_copy_last_pos_kernel")
		prof.end("cplist_copy_last_pos_kernel")

	def speak(self):
		print("CuPy cluster pair list: npa=%d cluster=%d width=%d rebuilds=%d checks=%d "
		      "cells=%d clusters_per_cell=%d exclusions=%s excl_width=%d"
		      % (int(self.info.npa), int(self.cluster_size), int(self.nmax),
		         int(self.rebuild_count), int(self.check_count), int(self.clist.ncell),
		         int(self.clusters_per_cell), ",".join(self.exclusion or ()),
		         int(self.excl_width)))


# Public class name used by pygamd.plist.neighbor.
cluster_pair = cplist
