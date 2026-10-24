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
# CuPy-native cell list (MVP).  Pure cupy.RawKernel; NO numba.  This is the main
# line: pygamd.plist.neighbor -> CuPy nlist -> CuPy clist.  The legacy Numba
# implementation has been moved to archive_legacy/pygamd/plists/clist.py.
#
# Scope (MVP): orthorhombic periodic box only.  Cells are sized to the search
# cutoff (rcut + rbuff), so all interaction partners of a particle live in the
# 27 cells adjacent to (and including) its own cell.  Cell capacity overflow is
# detected on device and the per-cell buffer is grown and the build retried.
#
# dtypes: particle index int32, cell size/offset/dim/adjacency int32, pos float32.
# cupy is imported lazily (kernels compiled on first build); importing this
# module never requires a CUDA driver.
# -----------------------------------------------------------------------------

import math
import numpy as np

from pygamd.backend.cupy import import_cupy
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace


# situation[0] = (nan particle index + 1), [1] = (out-of-box index + 1),
# [2] = max cell occupancy seen (for capacity overflow / growth).
_CLIST_SRC = r'''
extern "C" __global__ void clist_build(
		const int   npa,
		const float* __restrict__ pos,
		const int*   __restrict__ dim,
		const float* __restrict__ box_low,
		const float* __restrict__ inv_width,
		const int    ncell,
		const int    nmax,
		int*         __restrict__ cell_size,
		int*         __restrict__ cell_members,
		int*         __restrict__ situation)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float pix = pos[i * 4 + 0];
		float piy = pos[i * 4 + 1];
		float piz = pos[i * 4 + 2];
		if (isnan(pix) || isnan(piy) || isnan(piz)) { situation[0] = i + 1; return; }

		float lx = box_low[0], ly = box_low[1], lz = box_low[2];
		if (pix < lx || pix >= -lx || piy < ly || piy >= -ly || piz < lz || piz >= -lz)
		{ situation[1] = i + 1; return; }

		int ix = (int)((pix - lx) * inv_width[0]);
		int iy = (int)((piy - ly) * inv_width[1]);
		int iz = (int)((piz - lz) * inv_width[2]);
		if (ix == dim[0]) ix = 0;
		if (iy == dim[1]) iy = 0;
		if (iz == dim[2]) iz = 0;

		int cell_id = iz + dim[2] * (iy + ix * dim[1]);
		if (cell_id < 0 || cell_id >= ncell) { situation[1] = i + 1; return; }

		int slot = atomicAdd(&cell_size[cell_id], 1);
		if (slot < nmax) cell_members[cell_id * nmax + slot] = i;
		else             atomicMax(&situation[2], slot + 1);
	}
}
'''


class clist:
	def __init__(self, info, rcut):
		# rcut here is the SEARCH cutoff (rcut + rbuff) passed by nlist.
		self.rcut = rcut
		self.info = info
		self.block_size = 64
		self.cp = None
		self._k_build = None

		b = info.box
		self.box = np.asarray([b[0], b[1], b[2]], dtype=np.float32)
		self._compute_geometry()

		# per-cell capacity (grown on overflow); start from average occupancy.
		nmax = self.info.npa // max(int(self.ncell), 1)
		self.nmax = int(nmax + 8 - (nmax & 7))
		self._alloc_device()

	# --- host geometry --------------------------------------------------------
	def _compute_geometry(self):
		self.dim = np.asarray([int(self.box[0] / self.rcut),
		                       int(self.box[1] / self.rcut),
		                       int(self.box[2] / self.rcut)], dtype=np.int32)
		if self.dim[0] < 3 or self.dim[1] < 3 or self.dim[2] < 3:
			raise RuntimeError(
				"CuPy-native MVP clist requires an orthorhombic box of at least "
				"3*(rcut+rbuff) per dimension (got box=%s, cell dim=%s)."
				% (self.box.tolist(), self.dim.tolist()))
		self.width = np.asarray([self.box[0] / float(self.dim[0]),
		                         self.box[1] / float(self.dim[1]),
		                         self.box[2] / float(self.dim[2])], dtype=np.float32)
		self.inv_width = (1.0 / self.width).astype(np.float32)
		self.box_low_boundary = (-self.box / 2.0).astype(np.float32)
		self.ncell = int(self.dim[0] * self.dim[1] * self.dim[2])
		self.cell_adj = self._build_cell_map()   # (ncell, 27) int32, -1 padded

	def _build_cell_map(self):
		dimx, dimy, dimz = int(self.dim[0]), int(self.dim[1]), int(self.dim[2])
		adj = np.full((self.ncell, 27), -1, dtype=np.int32)
		for x in range(dimx):
			for y in range(dimy):
				for z in range(dimz):
					cell_id = z + dimz * (y + x * dimy)
					neigh = []
					for nx in range(x - 1, x + 2):
						for ny in range(y - 1, y + 2):
							for nz in range(z - 1, z + 2):
								wx = nx % dimx
								wy = ny % dimy
								wz = nz % dimz
								neigh.append(wz + dimz * (wy + wx * dimy))
					# unique neighbor cells (handles small dims without double count)
					uniq = np.unique(np.asarray(neigh, dtype=np.int32))
					adj[cell_id, :uniq.shape[0]] = uniq
		return adj

	# --- device setup ---------------------------------------------------------
	def _ensure_cupy(self):
		if self.cp is None:
			self.cp = import_cupy("CuPy cell list")
		if self._k_build is None:
			mod = self.cp.RawModule(code=_CLIST_SRC)
			self._k_build = mod.get_function("clist_build")
		return self.cp

	def _alloc_device(self):
		cp = self._ensure_cupy()
		self.d_dim = cp.asarray(self.dim)
		self.d_inv_width = cp.asarray(self.inv_width)
		self.d_box_low_boundary = cp.asarray(self.box_low_boundary)
		self.d_cell_adj = cp.asarray(self.cell_adj)            # (ncell, 27) int32
		self.d_cell_size = cp.zeros(self.ncell, dtype=cp.int32)
		self.d_cell_members = cp.zeros(self.ncell * self.nmax, dtype=cp.int32)
		self.d_situation = cp.zeros(3, dtype=cp.int32)

	def _update_geometry_if_box_changed(self):
		b = self.info.box
		if b[0] == self.box[0] and b[1] == self.box[1] and b[2] == self.box[2]:
			return
		self.box = np.asarray([b[0], b[1], b[2]], dtype=np.float32)
		self._compute_geometry()
		cp = self.cp
		self.d_dim = cp.asarray(self.dim)
		self.d_inv_width = cp.asarray(self.inv_width)
		self.d_box_low_boundary = cp.asarray(self.box_low_boundary)
		self.d_cell_adj = cp.asarray(self.cell_adj)
		self.d_cell_size = cp.zeros(self.ncell, dtype=cp.int32)
		self.d_cell_members = cp.zeros(self.ncell * self.nmax, dtype=cp.int32)

	# --- pos from the CuInfo CuPy device-data layer ---------------------------
	def _dev_pos(self):
		cu = getattr(self.info, "cu", None)
		if cu is not None:
			return cu.dev("pos")
		from pygamd.backend.cupy import as_cupy_array
		return as_cupy_array(self.info.d_pos, ravel=True)

	# --- build ----------------------------------------------------------------
	def calculate(self):
		prof.begin("clist_calculate_total")
		prof.count("clist_calculate_total")
		if _sp.ENABLED: _sp.inc("clist_build_count")
		cp = self._ensure_cupy()
		self._update_geometry_if_box_changed()
		pos = self._dev_pos()
		npa = int(self.info.npa)
		grid = (math.ceil(npa / self.block_size),)
		block = (self.block_size,)

		building = True
		while building:
			prof.begin("clist_cell_build_kernel")
			self.d_cell_size.fill(0)
			self.d_situation.fill(0)
			if _lt.LAUNCH: _lt.launch("clist_build_kernel")
			if _t.ENABLED: _t.gpu_begin("clist_build_kernel")
			self._k_build(grid, block,
			              (np.int32(npa), pos, self.d_dim, self.d_box_low_boundary,
			               self.d_inv_width, np.int32(self.ncell), np.int32(self.nmax),
			               self.d_cell_size, self.d_cell_members, self.d_situation))
			if _t.ENABLED: _t.gpu_end("clist_build_kernel")
			prof.end("clist_cell_build_kernel")
			prof.begin("clist_situation_sync")
			if _t.ENABLED: _t.host_begin("clist_status_copy")
			sit = cp.asnumpy(self.d_situation)
			if _t.ENABLED: _t.host_end("clist_status_copy")
			if _sp.ENABLED:
				_sp.inc("host_sync_count")            # device->host clist situation sync
				_sp.inc("host_sync_reason_clist")
			prof.end("clist_situation_sync")
			if int(sit[0]) > 0:
				raise RuntimeError('Error position NaN of particle ' + str(int(sit[0]) - 1) + ' !')
			if int(sit[1]) > 0:
				raise RuntimeError('Error particle ' + str(int(sit[1]) - 1) + ' out of box !')
			if int(sit[2]) > self.nmax:
				prof.count("clist_expand_count")
				nmax = int(sit[2])
				self.nmax = int(nmax + 8 - (nmax & 7))
				if _t.ENABLED: _t.host_begin("clist_realloc")
				self.d_cell_members = cp.zeros(self.ncell * self.nmax, dtype=cp.int32)
				if _t.ENABLED: _t.host_end("clist_realloc")
			else:
				building = False
		prof.end("clist_calculate_total")

	# --- fused-status opt-in fused-status helpers -------------------------------
	# Used ONLY by nlist's PYGAMD_NLIST_FUSED_STATUS path.  build_into() runs the
	# cell-build kernel ONCE into a caller-owned status buffer with NO host sync
	# and NO overflow re-loop; resolve_situation_host() applies the SAME safety
	# checks (NaN / out-of-box / cell overflow) from a host copy the caller made.
	# The default clist.calculate() above is untouched.
	def build_into(self, situation_slice):
		'''Run the cell-build kernel once into `situation_slice` (a (3,) device
		view), no host sync.  Caller must zero the slice beforehand and inspect it
		later via resolve_situation_host().'''
		cp = self._ensure_cupy()
		self._update_geometry_if_box_changed()
		pos = self._dev_pos()
		npa = int(self.info.npa)
		grid = (math.ceil(npa / self.block_size),)
		block = (self.block_size,)
		if _sp.ENABLED: _sp.inc("clist_build_count")
		prof.begin("clist_cell_build_kernel")
		self.d_cell_size.fill(0)
		if _lt.LAUNCH: _lt.launch("clist_build_kernel")
		if _t.ENABLED: _t.gpu_begin("clist_build_kernel")
		self._k_build(grid, block,
		              (np.int32(npa), pos, self.d_dim, self.d_box_low_boundary,
		               self.d_inv_width, np.int32(self.ncell), np.int32(self.nmax),
		               self.d_cell_size, self.d_cell_members, situation_slice))
		if _t.ENABLED: _t.gpu_end("clist_build_kernel")
		prof.end("clist_cell_build_kernel")

	def resolve_situation_host(self, sit):
		'''Given a host copy [nan, oob, maxocc] of a situation buffer: raise on
		NaN / out-of-box (same safety as calculate()); on cell overflow, grow the
		per-cell capacity and return True (caller must rebuild); else return False.
		Keeps every clist safety check -- nothing is skipped.'''
		if int(sit[0]) > 0:
			raise RuntimeError('Error position NaN of particle ' + str(int(sit[0]) - 1) + ' !')
		if int(sit[1]) > 0:
			raise RuntimeError('Error particle ' + str(int(sit[1]) - 1) + ' out of box !')
		if int(sit[2]) > self.nmax:
			prof.count("clist_expand_count")
			nmax = int(sit[2])
			self.nmax = int(nmax + 8 - (nmax & 7))
			if _t.ENABLED: _t.host_begin("clist_realloc")
			self.d_cell_members = self.cp.zeros(self.ncell * self.nmax, dtype=self.cp.int32)
			if _t.ENABLED: _t.host_end("clist_realloc")
			return True
		return False

	def speak(self):
		print("CuPy clist: ncell=%d nmax=%d dim=%s" % (self.ncell, self.nmax, self.dim.tolist()))
