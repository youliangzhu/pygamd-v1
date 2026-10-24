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
'''
CuPy-native DPD force (Groot-Warren dissipative particle dynamics).

Ported from GALAMOST DPDForce.cu/.cc + extern/saruprngCUDA.h.  One thread per
particle i, traverses the CuPy neighbor list, accumulates conservative +
dissipative + random pairwise forces into d_force / d_virial_potential.  The
virial and potential accounting intentionally follows pygamd-numba's DPD kernel
so A5B5-style scripts compare like-for-like.  cupy.RawKernel; no numba.

Per type-pair params (GALAMOST ToReal4):  [alpha, sigma, rcut^2, 1/rcut]
Per pair (rsq < rcut^2):
    rinv  = rsqrt(rsq);  omega = rinv - 1/rcut
    pair_eng = 0.25*alpha*omega^2*rsq*rcut
    gamma = 0.5*sigma^2
    Fc =  alpha*omega                                  (conservative)
    Fd = -dot(dx, dv) * gamma * omega^2                (dissipative)
    Fr =  sqrt(3) * sqrt(T) * sigma * omega * xi       (random; xi in [-1,1])
    Forcefac = Fc + Fd + Fr * (1/sqrt(dt))
    F_i += dx * Forcefac
Random number xi is drawn from Saru(seed, min(i,j), max(i,j)) with
seed = base_seed + timestep, matching pygamd-numba's index-based DPD RNG.
T_eff = sqrt(T) and rsq_dt = 1/sqrt(dt).

DPDGWVV support: RawDPD can read either the ordinary velocity field ('vel') or
the stable predictor velocity field ('velo') written by integration.dpdgwvv.
The public DPDGWVV entry sets velocity_source='velo'; ordinary DPD keeps 'vel'.
setConsiderDiameter remains intentionally unsupported in this lightweight branch.
'''

import math
import os
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.backend.cupy import import_cupy
from pygamd.backend.cupy.rng import SARU_SRC
from pygamd.backend.cupy.runtime import (
	VELOCITY_SOURCE_VEL, VELOCITY_SOURCE_VELO,
	normalize_dpd_velocity_source, require_cuinfo,
)


def _dpd_graph_dynamic_seed_enabled():
	v = os.environ.get("PYGAMD_DPD_GRAPH_DYNAMIC_SEED", "0").strip().lower()
	return v in ("1", "true", "on", "yes", "safe")


_DPD_SRC = SARU_SRC + r'''
extern "C" __global__ void dpd_force(
		const int    npa,
		const float* __restrict__ pos,
		const float* __restrict__ vel,
		const int*   __restrict__ tag,
		const float* __restrict__ params,
		const float* __restrict__ box,
		const int    ntypes,
		const int    nparams,
		const int*   __restrict__ neighbor_size,
		const int*   __restrict__ neighbor_list,
		const int    nlist_width,
		const unsigned int seed,
		const float  T_eff,
		const float  rsq_dt,
		float* __restrict__ force,
		float* __restrict__ virial_potential)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float pix = pos[i*4+0], piy = pos[i*4+1], piz = pos[i*4+2];
		int   ti  = (int)pos[i*4+3];
		float vix = vel[i*4+0], viy = vel[i*4+1], viz = vel[i*4+2];
		float bx = box[0], by = box[1], bz = box[2];
		float bxh = bx*0.5f, byh = by*0.5f, bzh = bz*0.5f;

		float f0 = force[i*3+0], f1 = force[i*3+1], f2 = force[i*3+2];
		float vir = virial_potential[i*2+0];
		float pot = virial_potential[i*2+1];

		int nsize = neighbor_size[i];
		for (int n = 0; n < nsize; n++)
		{
			int j = neighbor_list[i*nlist_width + n];
			float dx = pix - pos[j*4+0];
			float dy = piy - pos[j*4+1];
			float dz = piz - pos[j*4+2];
			if      (dx >=  bxh) dx -= bx; else if (dx < -bxh) dx += bx;
			if      (dy >=  byh) dy -= by; else if (dy < -byh) dy += by;
			if      (dz >=  bzh) dz -= bz; else if (dz < -bzh) dz += bz;
			float rsq = dx*dx + dy*dy + dz*dz;
			int tj = (int)pos[j*4+3];
			const float* pm = params + (ti + tj*ntypes) * nparams;
			float rcutsq = pm[2];
			if (rsq < rcutsq && rsq > 0.0f)
			{
				float dvx = vix - vel[j*4+0];
				float dvy = viy - vel[j*4+1];
				float dvz = viz - vel[j*4+2];
				float dot = dx*dvx + dy*dvy + dz*dvz;
				float alpha = pm[0], sigma = pm[1], rcutINV = pm[3];
				float rinv = rsqrtf(rsq);
				float omega = rinv - rcutINV;
				float pair_eng = 0.25f*alpha*omega*omega*rsq/rcutINV;
				float gamma = 0.5f*sigma*sigma;
				float Fc = alpha*omega;
				float Fd = -dot*gamma*omega*omega;
				// Match pygamd-numba: pair-symmetric Saru seeded by particle indices.
				unsigned int s2 = (i < j) ? (unsigned int)i : (unsigned int)j;
				unsigned int s3 = (i < j) ? (unsigned int)j : (unsigned int)i;
				Saru rng; saru3(&rng, seed, s2, s3);
				float xi = saru_f(&rng, -1.0f, 1.0f);
				float Fr = 1.7320508f * T_eff * sigma * omega * xi;
				float Forcefac = Fc + Fd + Fr*rsq_dt;
				f0 += dx*Forcefac; f1 += dy*Forcefac; f2 += dz*Forcefac;
				vir += (1.0f/6.0f)*rsq*Forcefac;
				pot += pair_eng;
			}
		}
		force[i*3+0] = f0; force[i*3+1] = f1; force[i*3+2] = f2;
		virial_potential[i*2+0] = vir;
		virial_potential[i*2+1] = pot;
	}
}
'''

# Graph-safe stochastic DPD variant.  It is intentionally identical to dpd_force
# except that seed and temperature are read from one-element device scalars.
# CUDA Graph records pointers, not values, so graph_prepare_replay() can update
# d_graph_seed and d_graph_temp_eff before each replay and preserve both
# Saru(base_seed+timestep, min_tag, max_tag) and VariantLinear temperature.
_DPD_SEEDPTR_SRC = (_DPD_SRC
	.replace('extern "C" __global__ void dpd_force(',
	         'extern "C" __global__ void dpd_force_seedptr(')
	.replace('const unsigned int seed,',
	         'const unsigned int* __restrict__ seed_ptr,')
	.replace('const float  T_eff,',
	         'const float* __restrict__ T_eff_ptr,')
	.replace('float* __restrict__ virial_potential)\n{',
	         'float* __restrict__ virial_potential)\n{\n\tconst unsigned int seed = seed_ptr[0];\n\tconst float T_eff = T_eff_ptr[0];', 1))


class dpd:
	'''CuPy-native DPD force backend.'''

	def __init__(self, info, nlist, rcut=1.0, temperature=1.0, seed=12345):
		self.info = info
		self.cu = require_cuinfo(info, "RawDPD force")
		self.nlist = nlist
		self.rcut = float(rcut)
		self.temp = temperature
		self.seed = int(seed)
		self.ntypes = int(info.ntypes)
		self.nparams = 4
		self.params = [[0.0] * self.nparams for _ in range(self.ntypes * self.ntypes)]
		self.params_changed = True
		self.block_size = 128
		self.cp = None
		self.kernel = None
		self.graph_seed_kernel = None
		self.d_params = None
		self.d_graph_seed = None
		self.d_graph_temp_eff = None
		# velocity source for dissipative DPD force.  Ordinary DPD uses d_vel;
		# integration.dpdgwvv switches this to d_velo before the force is captured.
		self.velocity_source = VELOCITY_SOURCE_VEL
		self._consider_diameter = False
		self._core_repulsion = 0.0
		# fused zero-force is NOT supported (DPD is additive over nlist)
		self.can_fuse_zero_force = False

	# --- params / options -----------------------------------------------------
	def setParams(self, type_i, type_j, alpha, sigma, r_cut=None, check=None):
		rc = float(self.rcut if r_cut is None else r_cut)
		ti = self.info.convert_name_to_id(type_i)
		tj = self.info.convert_name_to_id(type_j)
		row = [float(alpha), float(sigma), rc * rc, 1.0 / rc]
		idx1 = ti * self.ntypes + tj
		idx2 = tj * self.ntypes + ti
		self.params[idx1] = list(row)
		self.params[idx2] = list(row)
		if check is not None:
			check[idx1] = True
			check[idx2] = True
		self.params_changed = True

	def setT(self, temperature):
		self.temp = temperature

	def setVelocitySource(self, source):
		"""Select the velocity field used by the dissipative DPD force.

		This is an internal backend hook.  Public scripts should choose ordinary
		DPD by using force.dpd directly, or DPDGWVV by passing that force to
		integration.dpdgwvv(...).
		"""
		source = normalize_dpd_velocity_source(source)
		if source == VELOCITY_SOURCE_VELO and not self.cu.has_field(VELOCITY_SOURCE_VELO):
			raise RuntimeError("DPDGWVV velocity_source='velo' requires CuInfo.d_velo")
		self.velocity_source = source

	def _velocity_dev(self):
		return self.cu.dev(self.velocity_source)

	def setConsiderDiameter(self, con_dia, core_repulsion=0.0):
		self._core_repulsion = float(core_repulsion)
		if con_dia:
			raise NotImplementedError(
				"DPD setConsiderDiameter is partial: parameters are accepted but the "
				"diameter-correction kernel is not implemented this round.  Run DPD "
				"without diameter correction for now.")
		self._consider_diameter = False

	def _get_T(self, timestep):
		if hasattr(self.temp, "getValue"):
			return float(self.temp.getValue(timestep))
		return float(self.temp)

	# --- compute ---------------------------------------------------------------
	def set_fuse_zero_force(self, flag):
		# DPD never overwrites force; it is always additive.
		pass

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("RawDPD force")
			self.kernel = self.cp.RawKernel(_DPD_SRC, "dpd_force")
			self.graph_seed_kernel = self.cp.RawKernel(_DPD_SEEDPTR_SRC, "dpd_force_seedptr")
		return self.cp

	def calculate(self, timestep):
		cp = self._ensure()
		prof.begin("rawdpd_total")
		prof.count("rawdpd_total")
		if _sp.ENABLED:
			_sp.inc("dpd_compute_count")
			# DPD always computes conservative+dissipative+random in one kernel (Saru RNG)
			_sp.note("RawDPD: conservative+dissipative+random in one kernel per step; "
			         "Saru pair RNG keyed by (seed+timestep, min/max tag); "
			         "velocity_source=%s" % self.velocity_source)
		if self.params_changed:
			self.d_params = cp.asarray(np.asarray(self.params, dtype=np.float32).ravel())
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False
		self.nlist.calculate(timestep)
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		nlist_width = int(d_nlist.shape[1])
		pos = self.cu.dev("pos")
		vel = self._velocity_dev()
		tag = self.cu.dev("tag")
		box = self.cu.dev("box")
		force = self.cu.dev("force")
		vp = self.cu.dev("virial_potential")
		dt = float(self.info.dt)
		rsq_dt = (1.0 / math.sqrt(dt)) if dt > 1.0e-7 else 0.0
		T_eff = math.sqrt(max(self._get_T(timestep), 0.0))
		seed = np.uint32((self.seed + int(timestep)) & 0xffffffff)
		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		if _lt.LAUNCH: _lt.launch("dpd_force")
		if _t.ENABLED: _t.gpu_begin("dpd_force")
		self.kernel(grid, block,
		            (np.int32(self.info.npa), pos, vel, tag, self.d_params, box,
		             np.int32(self.ntypes), np.int32(self.nparams),
		             d_nsize, d_nlist.ravel(), np.int32(nlist_width),
		             seed, np.float32(T_eff), np.float32(rsq_dt), force, vp))
		if _t.ENABLED: _t.gpu_end("dpd_force")
		prof.end("rawdpd_total")

	# --- CUDA Graph/2L CUDA-graph helpers --------------------------------------
	# The graph path allowed conservative DPD (sigma == 0) because seed is irrelevant.
	# The stochastic graph path adds an opt-in graph-safe stochastic DPD path.  It preserves the current
	# Saru pair RNG semantics by reading base_seed+timestep from a device scalar that
	# is updated before every graph replay.  Default eager calculate() is unchanged.
	graph_name = "dpd"

	def _graph_sigma_zero(self):
		try:
			for row in self.params:
				if abs(float(row[1])) > 1.0e-12:
					return False
			return True
		except Exception:
			return False

	def _graph_temp_is_variant(self):
		return hasattr(self.temp, "getValue")

	def _graph_dynamic_seed_allowed(self):
		return _dpd_graph_dynamic_seed_enabled()

	def _graph_uses_dynamic_seed(self):
		return (not self._graph_sigma_zero()) and self._graph_dynamic_seed_allowed()

	def _ensure_graph_seed(self):
		cp = self._ensure()
		if self.d_graph_seed is None:
			self.d_graph_seed = cp.zeros((1,), dtype=cp.uint32)
		return self.d_graph_seed

	def _ensure_graph_temp_eff(self):
		cp = self._ensure()
		if self.d_graph_temp_eff is None:
			self.d_graph_temp_eff = cp.zeros((1,), dtype=cp.float32)
		return self.d_graph_temp_eff

	def graph_ready(self):
		base_ready = (self.cp is not None and self.kernel is not None and self.d_params is not None
		              and not self.params_changed)
		if not base_ready:
			return False
		if self._graph_sigma_zero():
			return True
		if not self._graph_dynamic_seed_allowed():
			return False
		self._ensure_graph_seed()
		self._ensure_graph_temp_eff()
		return (self.graph_seed_kernel is not None and self.d_graph_seed is not None
		        and self.d_graph_temp_eff is not None)

	def graph_unready_reason(self):
		if not self._graph_sigma_zero():
			if not self._graph_dynamic_seed_allowed():
				return "dpd_dynamic_seed"
			if self.graph_seed_kernel is None:
				return "dpd_graph_seed_kernel_missing"
		if self.velocity_source == VELOCITY_SOURCE_VELO and not self.cu.has_field(VELOCITY_SOURCE_VELO):
			return "dpd_velo_missing"
		if self.cp is None or self.kernel is None or self.d_params is None:
			return "dpd_not_initialized"
		if self.params_changed:
			return "dpd_params_changed"
		return "unsupported"

	def graph_prepare_replay(self, timestep):
		"""Update graph-dynamic scalar inputs immediately before graph.launch().

		For stochastic DPD this preserves the eager RNG key exactly:
		    seed = base_seed + timestep
		The graph captures device pointers to d_graph_seed and d_graph_temp_eff; replay
		reads the current values.  This is opt-in via PYGAMD_DPD_GRAPH_DYNAMIC_SEED=1
		for stochastic DPD and does not affect the default eager path.
		"""
		if not self._graph_uses_dynamic_seed():
			return
		d_seed = self._ensure_graph_seed()
		d_temp_eff = self._ensure_graph_temp_eff()
		seed = np.uint32((int(self.seed) + int(timestep)) & 0xffffffff)
		T_eff = np.float32(math.sqrt(max(float(self._get_T(timestep)), 0.0)))
		# CuPy fill enqueues tiny device operations/copies but avoids freezing by-value
		# seed/temperature scalars in the captured graph.  Future work can replace this
		# with graph node parameter updates or in-graph device counters.
		d_seed.fill(seed)
		d_temp_eff.fill(T_eff)
		if _sp.ENABLED:
			_sp.inc("cuda_graph_dpd_dynamic_seed_update_count")
			_sp.inc("cuda_graph_temperature_scalar_update_count")
			_sp.inc("cuda_graph_dpd_temperature_scalar_update_count")
			if self._graph_temp_is_variant():
				_sp.inc("cuda_graph_variant_temperature_update_count")
				_sp.inc("cuda_graph_dpd_variant_temperature_update_count")

	def graph_signature(self):
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		mode = "dpd_dynamic_seed" if self._graph_uses_dynamic_seed() else "dpd_conservative"
		seed_ptr_id = id(self.d_graph_seed) if self._graph_uses_dynamic_seed() else 0
		temp_ptr_id = id(self.d_graph_temp_eff) if self._graph_uses_dynamic_seed() else 0
		try:
			T0 = float(self._get_T(0))
			T_sig = ("variant_temperature", id(self.temp), T0) if self._graph_temp_is_variant() else ("static_temperature", T0)
		except Exception:
			T_sig = "unknown"
		return (mode, self.velocity_source, id(self.info.cu.dev("pos")), id(self._velocity_dev()),
		        id(self.info.cu.dev("tag")), id(self.info.cu.dev("force")),
		        id(self.info.cu.dev("virial_potential")), id(self.d_params),
		        id(d_nlist), id(d_nsize), int(d_nlist.shape[1]), int(self.ntypes),
		        int(self.nparams), int(self.info.npa), float(self.info.dt), int(self.seed),
		        T_sig, seed_ptr_id, temp_ptr_id)

	def graph_force_add(self):
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		nlist_width = int(d_nlist.shape[1])
		pos = self.cu.dev("pos")
		vel = self._velocity_dev()
		tag = self.cu.dev("tag")
		box = self.cu.dev("box")
		force = self.cu.dev("force")
		vp = self.cu.dev("virial_potential")
		dt = float(self.info.dt)
		rsq_dt = (1.0 / math.sqrt(dt)) if dt > 1.0e-7 else 0.0
		grid = (int(self.nblocks),)
		block = (int(self.block_size),)

		if self._graph_uses_dynamic_seed():
			self._ensure_graph_seed()
			self._ensure_graph_temp_eff()
			self.graph_seed_kernel(grid, block,
			            (np.int32(self.info.npa), pos, vel, tag, self.d_params, box,
			             np.int32(self.ntypes), np.int32(self.nparams),
			             d_nsize, d_nlist.ravel(), np.int32(nlist_width),
			             self.d_graph_seed, self.d_graph_temp_eff, np.float32(rsq_dt), force, vp))
			if _sp.ENABLED:
				_sp.inc("cuda_graph_dpd_dynamic_seed_force_count")
			return

		# Conservative-only graph path: sigma == 0 makes seed/T_eff irrelevant, but the
		# original kernel ABI still requires scalar values.
		seed = np.uint32(int(self.seed) & 0xffffffff)
		self.kernel(grid, block,
		            (np.int32(self.info.npa), pos, vel, tag, self.d_params, box,
		             np.int32(self.ntypes), np.int32(self.nparams),
		             d_nsize, d_nlist.ravel(), np.int32(nlist_width),
		             seed, np.float32(0.0), np.float32(rsq_dt), force, vp))
