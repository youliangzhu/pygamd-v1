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
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with PYGAMD are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''

# CuPy-native particle sort (tinker).  This round implements a DETERMINISTIC
# GPU CELL-KEY sort (particles reordered by their cell index in row-major cell
# order) -- a partial of the full GALAMOST Hilbert-curve traversal, which is a
# follow-up.  The whole hot sort path is now device-side:
#   - cell keys are computed from CuInfo.d_pos with CuPy,
#   - argsort uses a deterministic composite key on GPU,
#   - CuInfo.reorder() rebuilds d_rtag and topology indices on GPU,
#   - every particle_set's device member list is rebuilt from d_tag on GPU,
#   - every neighbor list is forced to rebuild on the next step.
# Tag-based topology (and dump.xml topology output) is sort-invariant.

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace


class sort_order:
	'''Periodic spatial sort of particles (cell-key; Hilbert traversal is a
	follow-up).  app.add not required -- dynamics(sort=True) drives it.'''

	def __init__(self, info):
		self.info = info
		self.cu = getattr(info, "cu", None)
		if self.cu is None:
			raise RuntimeError("Sort requires a CuInfo snapshot (backend='cupy').")
		self.grid = None          # optional (mx,my,mz); default derived from box
		self.sort_count = 0

	def setGrid(self, mx, my, mz):
		self.grid = (int(mx), int(my), int(mz))

	def _cell_dims(self):
		if self.grid is not None:
			return self.grid
		b = self.info.box
		# default: ~unit cells, capped, >=1 per axis (deterministic, box-derived)
		return (max(1, min(int(b[0]), 256)),
		        max(1, min(int(b[1]), 256)),
		        max(1, min(int(b[2]), 256)))

	def _cell_key_device(self):
		'''Return deterministic row-major cell keys as a CuPy device array.

		This replaces the previous host-key path (`cu.to_numpy("pos")` +
		NumPy cell-key calculation).  Positions are already kept in CuInfo.d_pos,
		so the key calculation can stay entirely on the current CUDA device.
		'''
		cp = self.cu.cp
		pos = self.cu.d_pos[:, :3]
		b = cp.asarray(self.info.box[:3], dtype=cp.float32)
		dim = cp.asarray(self._cell_dims(), dtype=cp.int32)
		low = -0.5 * b
		frac = (pos - low) / b
		idx = cp.floor(frac * dim).astype(cp.int32)
		idx = cp.clip(idx, 0, dim - 1)
		# row-major cell id: iz + dimz*(iy + ix*dimy)
		cell = idx[:, 2] + dim[2] * (idx[:, 1] + idx[:, 0] * dim[1])
		return cell.astype(cp.int64, copy=False)

	def _sort_order_device(self):
		'''Return a CuPy device permutation for the current particle order.

		CuPy's sort stability can vary by version/backend.  To preserve the old
		`np.argsort(..., kind="stable")` semantics, sort a composite key:
			cell_key * npa + current_index
		This orders by cell and keeps same-cell particles in current-index order.
		'''
		cp = self.cu.cp
		npa = int(self.info.npa)
		cell = self._cell_key_device()
		idx0 = cp.arange(npa, dtype=cp.int64)
		key = cell * cp.int64(npa) + idx0
		return cp.argsort(key).astype(cp.int32, copy=False)

	def calculate(self, timestep):
		prof.begin("tinker_sort_total")
		prof.count("tinker_sort_total")
		if _sp.ENABLED:
			_sp.inc("sort_count")
			_sp.note("tinker sort is GPU-key: CuPy cell-key + deterministic composite-key argsort; "
			         "CuInfo.reorder rebuilds rtag/topology on device; full Hilbert traversal deferred")
		if _lt.LAUNCH: _lt.launch("sort_launches_or_host_calls")
		if _t.ENABLED: _t.host_begin("sort")
		order = self._sort_order_device()
		if _t.ENABLED: _t.host_end("sort")
		# 1) reorder all per-particle CuPy arrays + rebuild rtag + refresh topology
		if _lt.LAUNCH: _lt.launch("topology_refresh_calls")
		if _t.ENABLED: _t.host_begin("topology_refresh")
		self.cu.reorder(order)
		# 2) rebuild every group's member list (tag-based)
		for ps in self.info.particle_set:
			ps.sort()
		# 3) force every neighbor list to rebuild on the next calculate()
		for nb in self.info.plist:
			data = getattr(nb, "data", None)
			if data is not None:
				data._built = False
		if _t.ENABLED: _t.host_end("topology_refresh")
		self.sort_count += 1
		prof.end("tinker_sort_total")
