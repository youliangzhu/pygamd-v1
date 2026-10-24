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

# -----------------------------------------------------------------------------
# RawKernel backend for custom non-bonded pair potentials (phase 2).
#
# This module wraps a user supplied CUDA __device__ pair potential (described by
# pygamd.potentials.RawPairPotential) into complete pair-force kernels using
# cupy.RawKernel.  It reuses the *existing* pygamd device arrays and neighbor
# list; it does NOT rebuild the neighbor list nor copy any core array to host.
#
# Backend optimizations (backend-only; user API unchanged):
#   * two kernels are built from the same user device function:
#       - FULL      : writes force + virial + potential (unchanged numerics),
#       - FORCE_ONLY: writes force only, never touches virial_potential, and
#                     does not accumulate the potential term (fp[1] is left to
#                     the compiler to dead-code-eliminate when the user function
#                     is __forceinline__).
#   * each timestep the backend picks FORCE_ONLY unless the current
#     info.compute_properties asks for potential / pressure / virial /
#     stress_tensor (then it conservatively uses FULL).
#   * optional --use_fast_math NVRTC option driven by potential.fast_math.
#   * cupy wrappers of the persistent device arrays are cached to cut per-step
#     Python overhead (re-wrapped only if the underlying array object changes).
#   * a very light opt-in profiler (env PYGAMD_RAWPAIR_PROFILE=1) reports how
#     often each kernel is taken, without changing physical results.
#
# Important: cupy is imported lazily, only when this backend is actually used.
# `import pygamd` never depends on cupy.
# -----------------------------------------------------------------------------

import os
import math
import numpy as np
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.forces import potentials as _potentials
from pygamd.backend.cupy import import_cupy, as_cupy_array


# Shared neighbor-traversal body.  %(extra_decl)s / %(extra_acc)s / %(extra_store)s
# splice in the virial/potential parts only for the FULL kernel.
_KERNEL_TEMPLATE = r'''
extern "C" __global__ void %(kname)s(
		const int    npa,
		const float* __restrict__ pos,
		const float* __restrict__ params,
		const float* __restrict__ box,
		const int    ntypes,
		const int    nparams,
		const int*   __restrict__ neighbor_size,
		const int*   __restrict__ neighbor_list,
		const int    nlist_width,
		float*       __restrict__ force%(extra_args)s)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float box0 = box[0];
		float box1 = box[1];
		float box2 = box[2];
		float box0_half = box0 * 0.5f;
		float box1_half = box1 * 0.5f;
		float box2_half = box2 * 0.5f;

		float pix = pos[i * 4 + 0];
		float piy = pos[i * 4 + 1];
		float piz = pos[i * 4 + 2];
		int   ti  = (int)pos[i * 4 + 3];

		// init in ADD mode reads the existing force (relies on a prior global
		// zero_force_vp); in WRITE/fused mode it starts from 0 and overwrites,
		// so a particle with no neighbors still ends up with force = 0.
		float result0 = %(f0)s;
		float result1 = %(f1)s;
		float result2 = %(f2)s;
%(extra_decl)s
		int nsize = neighbor_size[i];
		for (int ni = 0; ni < nsize; ni++)
		{
			int   j   = neighbor_list[i * nlist_width + ni];
			float dp0 = pix - pos[j * 4 + 0];
			float dp1 = piy - pos[j * 4 + 1];
			float dp2 = piz - pos[j * 4 + 2];
			int   tj  = (int)pos[j * 4 + 3];

			if      (dp0 >=  box0_half) dp0 -= box0;
			else if (dp0 <  -box0_half) dp0 += box0;

			if      (dp1 >=  box1_half) dp1 -= box1;
			else if (dp1 <  -box1_half) dp1 += box1;

			if      (dp2 >=  box2_half) dp2 -= box2;
			else if (dp2 <  -box2_half) dp2 += box2;

			float rsq = dp0 * dp0 + dp1 * dp1 + dp2 * dp2;
			const float* pms = params + (ti + tj * ntypes) * nparams;

			float fp[2];
			fp[0] = 0.0f;
			fp[1] = 0.0f;
			%(uname)s(rsq, pms, fp);

			result0 += fp[0] * dp0;
			result1 += fp[0] * dp1;
			result2 += fp[0] * dp2;
%(extra_acc)s
		}

		force[i * 3 + 0] = result0;
		force[i * 3 + 1] = result1;
		force[i * 3 + 2] = result2;
%(extra_store)s
	}
}
'''

# pieces that make the FULL kernel compute virial + potential
_FULL_EXTRA_ARGS = ",\n\t\tfloat* __restrict__ virial_potential,\n\t\tconst float one_sixth"
_FULL_EXTRA_ACC = (
	"\t\t\tresult3 += one_sixth * rsq * fp[0];\n"
	"\t\t\tresult4 += fp[1] * 0.5f;\n"
)
_FULL_EXTRA_STORE = (
	"\t\tvirial_potential[i * 2 + 0] = result3;\n"
	"\t\tvirial_potential[i * 2 + 1] = result4;\n"
)


def _kernel_pieces(full, write):
	'''Build the template substitution dict for one of the 4 kernel variants.

	full  : True -> compute/store virial + potential (FULL); False -> FORCE-only.
	write : True -> overwrite mode (fused zero, start from 0);
	        False -> additive mode (relies on a prior global zero_force_vp).
	'''
	if write:
		f0 = f1 = f2 = "0.0f"
		vp_decl = "\t\tfloat result3 = 0.0f;\n\t\tfloat result4 = 0.0f;\n"
	else:
		f0 = "force[i * 3 + 0]"
		f1 = "force[i * 3 + 1]"
		f2 = "force[i * 3 + 2]"
		vp_decl = (
			"\t\tfloat result3 = virial_potential[i * 2 + 0];\n"
			"\t\tfloat result4 = virial_potential[i * 2 + 1];\n"
		)
	if full:
		return {
			"f0": f0, "f1": f1, "f2": f2,
			"extra_args": _FULL_EXTRA_ARGS,
			"extra_decl": vp_decl,
			"extra_acc": _FULL_EXTRA_ACC,
			"extra_store": _FULL_EXTRA_STORE,
		}
	return {
		"f0": f0, "f1": f1, "f2": f2,
		"extra_args": "",
		"extra_decl": "",
		"extra_acc": "",
		"extra_store": "",
	}


def _import_cupy():
	'''Lazy import of cupy via the central CuPy bridge (clear error message).'''
	return import_cupy("RawKernel pair potential")


def _profile_enabled():
	v = os.environ.get("PYGAMD_RAWPAIR_PROFILE", "")
	return v not in ("", "0", "false", "False", "no", "off")


def _pair_charge_force_full_enabled():
	v = os.environ.get("PYGAMD_PAIR_CHARGE_FORCE_FULL", "")
	return v not in ("", "0", "false", "False", "no", "off", "NO", "OFF")


def _read_block_size(default):
	# Opt-in generic block-size override (audit/sweep only).
	# Applies uniformly to all RawPair potentials; default keeps prior behavior.
	v = os.environ.get("PYGAMD_RAWPAIR_BLOCK_SIZE", "")
	if v == "":
		return default
	try:
		b = int(v)
	except Exception:
		return default
	if b < 32 or b > 1024 or (b % 32) != 0:
		return default
	return b


_PAIR_CHARGE_KERNEL_TEMPLATE = r'''
%(device_code)s

extern "C" __global__ void %(kname)s(
		const int    npa,
		const float* __restrict__ pos,
		const float* __restrict__ charge,
		const float* __restrict__ params,
		const float* __restrict__ box,
		const int    ntypes,
		const int    nparams,
		const int*   __restrict__ neighbor_size,
		const int*   __restrict__ neighbor_list,
		const int    nlist_width,
		float*       __restrict__ force%(extra_args)s)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < npa)
	{
		float box0 = box[0];
		float box1 = box[1];
		float box2 = box[2];
		float box0_half = box0 * 0.5f;
		float box1_half = box1 * 0.5f;
		float box2_half = box2 * 0.5f;

		float pix = pos[i * 4 + 0];
		float piy = pos[i * 4 + 1];
		float piz = pos[i * 4 + 2];
		float qi = charge[i];
		int ti = (int)pos[i * 4 + 3];

		float result0 = %(f0)s;
		float result1 = %(f1)s;
		float result2 = %(f2)s;
%(extra_decl)s

		int nsize = neighbor_size[i];
		for (int ni = 0; ni < nsize; ni++)
		{
			int j = neighbor_list[i * nlist_width + ni];
			float dp0 = pix - pos[j * 4 + 0];
			float dp1 = piy - pos[j * 4 + 1];
			float dp2 = piz - pos[j * 4 + 2];
			float qj = charge[j];
			int tj = (int)pos[j * 4 + 3];

			if      (dp0 >=  box0_half) dp0 -= box0;
			else if (dp0 <  -box0_half) dp0 += box0;
			if      (dp1 >=  box1_half) dp1 -= box1;
			else if (dp1 <  -box1_half) dp1 += box1;
			if      (dp2 >=  box2_half) dp2 -= box2;
			else if (dp2 <  -box2_half) dp2 += box2;

			float rsq = dp0 * dp0 + dp1 * dp1 + dp2 * dp2;
			const float* pms = params + (ti + tj * ntypes) * nparams;

			float fp[2];
			fp[0] = 0.0f;
			fp[1] = 0.0f;
			%(uname)s(rsq, qi, qj, pms, fp);

			result0 += fp[0] * dp0;
			result1 += fp[0] * dp1;
			result2 += fp[0] * dp2;
%(extra_acc)s
		}

		force[i * 3 + 0] = result0;
		force[i * 3 + 1] = result1;
		force[i * 3 + 2] = result2;
%(extra_store)s
	}
}
'''


def _as_pair_charge_potential(func):
	if isinstance(func, str):
		if func in ("lj_coulomb", "lj_coulomb_erfc", "coulomb_erfc"):
			return _potentials.get_pair_charge_builtin(func)
		raise RuntimeError(
			"CuPy nonbonded_c currently supports built-in 'lj_coulomb', "
			"'lj_coulomb_erfc', 'coulomb_erfc' or "
			"pygamd.potentials.pair_charge(...); got %r." % (func,))
	if getattr(func, "kind", None) == "pair_charge":
		return func
	raise RuntimeError(
		"CuPy nonbonded_c requires a built-in charged pair string or a charged pair "
		"potential from pygamd.potentials.pair_charge(...).")


class pair_charge:
	'''CuPy RawKernel backend for charged non-bonded pair potentials.

	This mirrors pygamd-numba force.nonbonded_c.  Built-in func='lj_coulomb'
	uses user parameters [epsilon, sigma, alpha, epsilonr, rcut].
	For PPPM real space, func='lj_coulomb_erfc' uses
	[epsilon, sigma, alpha, epsilonr, rcut, kappa].
	Custom potentials are built with pygamd.potentials.pair_charge(...).
	'''

	graph_name = "pair_charge"

	def __init__(self, info, nlist, func, block_size=128):
		self.info = info
		self.cu = getattr(info, "cu", None)
		self.nlist = nlist
		self.block_size = _read_block_size(block_size)
		self.cp = None
		self.kernel_full_add = None
		self.kernel_force_add = None
		self.kernel_full_write = None
		self.kernel_force_write = None
		self.params_changed = True
		self._wrap_cache = {}
		self.can_fuse_zero_force = True
		self.fuse_zero_force_this_step = False
		self.force_full = _pair_charge_force_full_enabled()

		self.potential = _as_pair_charge_potential(func)
		self.func = self.potential.name
		self.user_nparams = int(self.potential.nparams)
		self.nparams = int(getattr(self.potential, "device_nparams", self.user_nparams))
		self.param_style = getattr(self.potential, "param_style", None)
		self.fast_math = bool(getattr(self.potential, "fast_math", False))
		self.device_code = self.potential.code
		self.func_name = self.potential.name

		self.params = [[] for _ in range(self.info.ntypes * self.info.ntypes)]

	def set_fuse_zero_force(self, flag):
		self.fuse_zero_force_this_step = bool(flag)

	def _compile(self, kname, pieces):
		if self.cp is None:
			self.cp = _import_cupy()
		pieces = dict(pieces)
		pieces.update({
			"device_code": self.device_code,
			"kname": kname,
			"uname": self.func_name,
		})
		source = _PAIR_CHARGE_KERNEL_TEMPLATE % {
			"device_code": pieces["device_code"],
			"kname": pieces["kname"],
			"uname": pieces["uname"],
			"f0": pieces["f0"],
			"f1": pieces["f1"],
			"f2": pieces["f2"],
			"extra_args": pieces["extra_args"],
			"extra_decl": pieces["extra_decl"],
			"extra_acc": pieces["extra_acc"],
			"extra_store": pieces["extra_store"],
		}
		options = ('--use_fast_math',) if self.fast_math else ()
		try:
			kernel = self.cp.RawKernel(source, kname, options=options)
			kernel.compile()
			return kernel
		except Exception as e:
			raise RuntimeError(
				"Failed to compile CuPy charged pair potential '%s'.\n"
				"Expected device signature:\n"
				"    __device__ void %s(float rsq, float qi, float qj, const float* param, float* fp)\n"
				"Original error: %s" % (self.func, self.func_name, repr(e))
			) from e

	def _build_one(self, suffix, full, write):
		kname = "pygamd_pair_charge_" + self.func_name + suffix
		return self._compile(kname, _kernel_pieces(full, write))

	def _build_kernels(self):
		self.kernel_full_add    = self._build_one("_full_add",    full=True,  write=False)
		self.kernel_force_add   = self._build_one("_force_add",   full=False, write=False)
		self.kernel_full_write  = self._build_one("_full_write",  full=True,  write=True)
		self.kernel_force_write = self._build_one("_force_write", full=False, write=True)

	def _wrap(self, key, arr):
		cached = self._wrap_cache.get(key)
		if cached is not None and cached[0] is arr:
			return cached[1]
		w = as_cupy_array(arr, ravel=True)
		self._wrap_cache[key] = (arr, w)
		return w

	def _dev(self, name, fallback_arr):
		if self.cu is not None:
			return self.cu.dev(name)
		return self._wrap(name, fallback_arr)

	def _need_virial_potential(self):
		if self.force_full:
			return True
		cp = getattr(self.info, "compute_properties", None)
		if not isinstance(cp, dict):
			return True
		return bool(
			cp.get("potential", False)
			or cp.get("pressure", False)
			or cp.get("virial", False)
			or cp.get("stress_tensor", False)
		)

	def setParams(self, type_i, type_j, param, check):
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id

		if self.param_style == "lj_coulomb":
			if len(param) != 5:
				raise RuntimeError("CuPy lj_coulomb expects 5 parameters [epsilon, sigma, alpha, epsilonr, rcut]; got %d." % len(param))
			epsilon = float(param[0])
			sigma = float(param[1])
			alpha = float(param[2])
			epsilonr = float(param[3])
			rcut = float(param[4])
			if epsilonr == 0.0:
				raise RuntimeError("CuPy lj_coulomb epsilonr must be non-zero.")
			p0 = 4.0 * epsilon * math.pow(sigma, 12)
			p1 = alpha * 4.0 * epsilon * math.pow(sigma, 6)
			stored = [p0, p1, 138.935 / epsilonr, rcut * rcut]
		elif self.param_style == "lj_coulomb_erfc":
			if len(param) != 6:
				raise RuntimeError("CuPy lj_coulomb_erfc expects 6 parameters [epsilon, sigma, alpha, epsilonr, rcut, kappa]; got %d." % len(param))
			epsilon = float(param[0])
			sigma = float(param[1])
			alpha = float(param[2])
			epsilonr = float(param[3])
			rcut = float(param[4])
			kappa = float(param[5])
			if epsilonr == 0.0:
				raise RuntimeError("CuPy lj_coulomb_erfc epsilonr must be non-zero.")
			if kappa <= 0.0:
				raise RuntimeError("CuPy lj_coulomb_erfc kappa must be positive.")
			p0 = 4.0 * epsilon * math.pow(sigma, 12)
			p1 = alpha * 4.0 * epsilon * math.pow(sigma, 6)
			stored = [p0, p1, 138.935 / epsilonr, kappa, rcut * rcut]
		elif self.param_style == "coulomb_erfc":
			if len(param) != 3:
				raise RuntimeError("CuPy coulomb_erfc expects 3 parameters [epsilonr, rcut, kappa]; got %d." % len(param))
			epsilonr = float(param[0])
			rcut = float(param[1])
			kappa = float(param[2])
			if epsilonr == 0.0:
				raise RuntimeError("CuPy coulomb_erfc epsilonr must be non-zero.")
			if kappa <= 0.0:
				raise RuntimeError("CuPy coulomb_erfc kappa must be positive.")
			stored = [138.935 / epsilonr, kappa, rcut * rcut]
		else:
			if len(param) != self.user_nparams:
				raise RuntimeError("CuPy charged pair potential '%s' expects %d parameters; got %d."
				                   % (self.potential.name, self.user_nparams, len(param)))
			stored = [float(x) for x in param]

		self.params[idx1] = list(stored)
		self.params[idx2] = list(stored)
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True

	def calculate(self, timestep):
		if getattr(self.info, "d_charge", None) is None:
			raise RuntimeError("CuPy nonbonded_c requires per-particle charges in the snapshot.")
		prof.begin("pair_charge_total")
		prof.count("pair_charge_total")
		if self.cp is None:
			self.cp = _import_cupy()
		if self.kernel_full_add is None:
			prof.begin("pair_charge_kernel_compile")
			self._build_kernels()
			prof.end("pair_charge_kernel_compile")
		if _sp.ENABLED:
			_sp.inc("pair_charge_compute_count")
			_sp.set_meta("pair_charge_block_size", int(self.block_size))

		if self.params_changed:
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = as_cupy_array(self.h_params, ravel=True)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False

		self.nlist.calculate(timestep)
		d_nlist = self.nlist.data.d_neighbor_list
		d_nlist_size = self.nlist.data.d_neighbor_size
		nlist_width = int(d_nlist.shape[1])

		pos = self._dev("pos", self.info.d_pos)
		charge = self._dev("charge", self.info.d_charge)
		box = self._dev("box", self.info.d_box)
		force = self._dev("force", self.info.d_force)

		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		common = (
			np.int32(self.info.npa),
			pos,
			charge,
			self.d_params,
			box,
			np.int32(self.info.ntypes),
			np.int32(self.nparams),
			d_nlist_size,
			d_nlist.ravel(),
			np.int32(nlist_width),
			force,
		)
		need_full = self._need_virial_potential()
		write = self.fuse_zero_force_this_step
		if _lt.LAUNCH:
			_lt.launch("pair_charge_force")
		if _t.ENABLED:
			_t.gpu_begin("pair_charge_force")
		if write:
			if need_full:
				virial = self._dev("virial_potential", self.info.d_virial_potential)
				self.kernel_full_write(grid, block, common + (virial, np.float32(1.0 / 6.0)))
			else:
				self.kernel_force_write(grid, block, common)
		else:
			if need_full:
				virial = self._dev("virial_potential", self.info.d_virial_potential)
				self.kernel_full_add(grid, block, common + (virial, np.float32(1.0 / 6.0)))
			else:
				self.kernel_force_add(grid, block, common)
		if _t.ENABLED:
			_t.gpu_end("pair_charge_force")
		prof.end("pair_charge_total")

	def graph_ready(self):
		return (self.cp is not None and self.kernel_force_write is not None
		        and self.kernel_force_add is not None and not self.params_changed)

	def graph_unready_reason(self):
		if self.cp is None or self.kernel_force_write is None or self.kernel_force_add is None:
			return "pair_charge_not_initialized"
		if self.params_changed:
			return "pair_charge_params_changed"
		return "unsupported"

	def graph_signature(self):
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		return (id(self.info.d_pos), id(self.info.d_charge), id(self.info.d_force),
		        id(self.d_params), id(d_nlist), id(d_nsize), int(d_nlist.shape[1]),
		        int(self.info.ntypes), int(self.info.npa))

	def _graph_common_force_args(self):
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		nlist_width = int(d_nlist.shape[1])
		pos = self._dev("pos", self.info.d_pos)
		charge = self._dev("charge", self.info.d_charge)
		box = self._dev("box", self.info.d_box)
		force = self._dev("force", self.info.d_force)
		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		common = (
			np.int32(self.info.npa),
			pos,
			charge,
			self.d_params,
			box,
			np.int32(self.info.ntypes),
			np.int32(self.nparams),
			d_nsize,
			d_nlist.ravel(),
			np.int32(nlist_width),
			force,
		)
		return grid, block, common

	def graph_force_write(self):
		grid, block, common = self._graph_common_force_args()
		self.kernel_force_write(grid, block, common)

	def graph_force_add(self):
		grid, block, common = self._graph_common_force_args()
		self.kernel_force_add(grid, block, common)


class pair:
	'''RawKernel backend mirroring forces.pair.pair for custom pair potentials.

	It exposes the same setParams(type_i, type_j, param, check) interface as the
	Numba pair backend so that forces.pair.pair can dispatch to it transparently.
	'''

	# how often (in calculate() calls) the profiler prints, when enabled
	_PROFILE_PERIOD = 1000

	def __init__(self, info, nlist, potential, block_size=128):
		if getattr(potential, "kind", None) != "pair_raw":
			raise RuntimeError(
				"Error: pygamd.force.nonbonded/pair requires a CuPy RawKernel pair "
				"potential built with pygamd.potentials.pair_raw / lj_raw / gem4_raw. "
				"The legacy Numba device-function pair path has been removed in the "
				"CuPy-native MVP. Example: "
				"func=pygamd.potentials.gem4_raw(fast_math=True).")
		self.info = info
		# CuPy-native device-data layer (MVP-1): when the snapshot provides a
		# CuInfo (info.cu), the core particle CuPy arrays come from it directly
		# instead of being re-wrapped from scattered Numba device arrays here.
		self.cu = getattr(info, "cu", None)
		self.nlist = nlist
		self.potential = potential
		# Default RawPair launch block size is 128 (GENERIC, applies to
		# every RawPair potential identically; NOT specialized by type/potential).
		# Validated via cross-engine data.log/RDF correctness (PASS) and trusted
		# performance: block=128 was clearly faster than 256 for GCM and LJ.
		# Override (and revert to 256) with PYGAMD_RAWPAIR_BLOCK_SIZE=256.
		self.block_size = _read_block_size(block_size)
		self.nparams = int(potential.nparams)
		self.fast_math = bool(getattr(potential, "fast_math", False))

		# same flatten layout as the Numba path: index = ti + tj*ntypes
		self.params = [[] for i in range(self.info.ntypes * self.info.ntypes)]
		self.params_changed = True

		# cupy + RawKernel are only touched the first time we actually compute.
		# Four kernel variants (Python-side selection, no per-thread branching):
		#   add   : additive  (relies on a prior global zero_force_vp)
		#   write : overwrite (fused zero-force, kernel writes the total directly)
		self.cp = None
		self.kernel_full_add = None
		self.kernel_force_add = None
		self.kernel_full_write = None
		self.kernel_force_write = None

		# per-step fused-zero flag, set by application via set_fuse_zero_force().
		# Default False -> additive (additive behavior preserved).
		self.can_fuse_zero_force = True
		self.fuse_zero_force_this_step = False

		# cache of cupy wrappers for persistent device arrays (id-guarded)
		self._wrap_cache = {}

		# lightweight profiler state (split by add/write x full/force)
		self._profile = _profile_enabled()
		self._calc_count = 0
		self._n_add_full = 0
		self._n_add_force = 0
		self._n_write_full = 0
		self._n_write_force = 0
		self._n_fused = 0

	# fused zero-force control: application toggles this each step
	def set_fuse_zero_force(self, flag):
		self.fuse_zero_force_this_step = bool(flag)

	# --- kernel compilation ---------------------------------------------------
	def _compile(self, kname, pieces):
		cp = self.cp
		source = self.potential.code + "\n" + (_KERNEL_TEMPLATE % pieces)
		options = ('--use_fast_math',) if self.fast_math else ()
		try:
			k = cp.RawKernel(source, kname, options=options)
			k.compile()  # force compilation now so errors surface immediately
			return k
		except Exception as e:
			raise RuntimeError(
				"Failed to compile RawKernel pair potential '%s' (kernel '%s', fast_math=%s).\n"
				"Please check that your CUDA device function signature is exactly:\n"
				"    __device__ void %s(float rsq, const float* param, float* fp)\n"
				"(  __device__ __forceinline__ void %s(...) is also supported  )\n"
				"Original error: %s"
				% (self.potential.name, kname, self.fast_math,
				   self.potential.name, self.potential.name, repr(e))
			) from e

	def _build_one(self, suffix, full, write):
		base = "pygamd_pair_" + self.potential.name
		kname = base + suffix
		pieces = _kernel_pieces(full, write)
		pieces["kname"] = kname
		pieces["uname"] = self.potential.name
		return self._compile(kname, pieces)

	def _build_kernels(self):
		# 4 variants, selected on the Python side (no runtime add/write branch)
		self.kernel_full_add    = self._build_one("_full_add",    full=True,  write=False)
		self.kernel_force_add   = self._build_one("_force_add",   full=False, write=False)
		self.kernel_full_write  = self._build_one("_full_write",  full=True,  write=True)
		self.kernel_force_write = self._build_one("_force_write", full=False, write=True)

	# --- parameter setting ----------------------------------------------------
	# NOTE: identical signature to forces.pair.pair.setParams so it can be
	#       dispatched to without changing the user-facing API.
	def setParams(self, type_i, type_j, param, check):
		if len(param) != self.nparams:
			raise RuntimeError(
				"Error, RawKernel pair potential '%s' expects %d parameters, got %d!"
				% (self.potential.name, self.nparams, len(param))
			)
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id
		self.params[idx1] = list(param)
		self.params[idx2] = list(param)
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True

	# --- kernel selection -----------------------------------------------------
	def _need_virial_potential(self):
		'''Conservatively decide whether this timestep needs virial/potential.

		Returns True (use FULL kernel) if any property that consumes
		d_virial_potential is requested for the current timestep.  Unknown or
		missing keys are treated as False; when in doubt we keep FULL.
		'''
		cp = getattr(self.info, "compute_properties", None)
		if not isinstance(cp, dict):
			return True  # cannot tell -> be safe
		return bool(
			cp.get("potential", False)
			or cp.get("pressure", False)
			or cp.get("virial", False)
			or cp.get("stress_tensor", False)
		)

	# --- cached cupy wrapper (zero-copy) --------------------------------------
	def _wrap(self, key, arr):
		cached = self._wrap_cache.get(key)
		if cached is not None and cached[0] is arr:
			return cached[1]
		w = as_cupy_array(arr, ravel=True)
		self._wrap_cache[key] = (arr, w)
		return w

	# --- CuInfo-preferred core-array accessor ---------------------------------
	def _dev(self, name, fallback_arr):
		'''Return a 1-D CuPy view of a core particle array.  Preferred path: ask
		the CuPy-native CuInfo device-data layer (centralized).  Transitional
		path (no CuInfo): zero-copy wrap of the Numba device array.  Both yield a
		CuPy array for the SAME RawKernel (the physics path never branches).'''
		if self.cu is not None:
			return self.cu.dev(name)
		return self._wrap(name, fallback_arr)

	# --- force computation ----------------------------------------------------
	def calculate(self, timestep):
		prof.begin("pair_total")
		prof.count("pair_total")
		if _sp.ENABLED:
			_sp.inc("pair_compute_count")
			_sp.set_meta("pair_block_size", int(self.block_size))
		# lazy: import cupy and compile kernels on first real use only
		if self.cp is None:
			self.cp = _import_cupy()
		if self.kernel_full_add is None:
			prof.begin("pair_kernel_compile")
			self._build_kernels()
			prof.end("pair_kernel_compile")
		cp = self.cp

		if self.params_changed:
			# flatten params: shape [ntypes*ntypes, nparams], same as Numba path
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = as_cupy_array(self.h_params, ravel=True)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False

		self.nlist.calculate(timestep)
		d_nlist = self.nlist.data.d_neighbor_list      # CuPy int32 [npa, width]
		d_nlist_size = self.nlist.data.d_neighbor_size # CuPy int32 [npa]
		nlist_width = int(d_nlist.shape[1])

		# core particle arrays via the CuInfo CuPy device-data layer (centralized).
		# pos/box/force/virial are owned by CuInfo; the neighbor-list arrays are
		# already CuPy (CuPy-native nlist) and are consumed directly -- NO Numba
		# wrap.  d_neighbor_list may be reallocated on overflow, so ravel each step.
		prof.begin("pair_wrap_arrays")
		prof.count("pair_wrap_arrays")
		pos    = self._dev("pos", self.info.d_pos)
		box    = self._dev("box", self.info.d_box)
		force  = self._dev("force", self.info.d_force)
		nlist  = d_nlist.ravel()
		nsize  = d_nlist_size
		prof.end("pair_wrap_arrays")

		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		common = (
			np.int32(self.info.npa),
			pos,
			self.d_params,
			box,
			np.int32(self.info.ntypes),
			np.int32(self.nparams),
			nsize,
			nlist,
			np.int32(nlist_width),
			force,
		)

		need_full = self._need_virial_potential()
		write = self.fuse_zero_force_this_step
		prof.count("pair_need_full_true" if need_full else "pair_need_full_false")
		prof.count("pair_fused_true" if write else "pair_fused_false")

		if _lt.LAUNCH: _lt.launch("pair_force")
		if _t.ENABLED: _t.gpu_begin("pair_force")
		if write:
			self._n_fused += 1
			if need_full:
				virial = self._dev("virial_potential", self.info.d_virial_potential)
				prof.begin("pair_full_write")
				self.kernel_full_write(grid, block, common + (virial, np.float32(1.0 / 6.0)))
				prof.end("pair_full_write")
				self._n_write_full += 1
			else:
				# force-only fused: kernel overwrites force; virial_potential is
				# intentionally left untouched (not needed this step, will be
				# overwritten by the next full_write step).
				prof.begin("pair_force_write")
				self.kernel_force_write(grid, block, common)
				prof.end("pair_force_write")
				self._n_write_force += 1
		else:
			if need_full:
				virial = self._dev("virial_potential", self.info.d_virial_potential)
				prof.begin("pair_full_add")
				self.kernel_full_add(grid, block, common + (virial, np.float32(1.0 / 6.0)))
				prof.end("pair_full_add")
				self._n_add_full += 1
			else:
				prof.begin("pair_force_add")
				self.kernel_force_add(grid, block, common)
				prof.end("pair_force_add")
				self._n_add_force += 1

		if _t.ENABLED: _t.gpu_end("pair_force")
		prof.end("pair_total")

		# --- opt-in lightweight profiler (no effect on physics) ---------------
		if self._profile:
			self._calc_count += 1
			if self._calc_count % self._PROFILE_PERIOD == 0:
				print("info : [pair profile] name=%s fast_math=%s calls=%d "
				      "add_full=%d add_force=%d write_full=%d write_force=%d fused=%d"
				      % (self.potential.name, self.fast_math, self._calc_count,
				         self._n_add_full, self._n_add_force,
				         self._n_write_full, self._n_write_force, self._n_fused))

	# --- CUDA-graph helpers (launch-only; no nlist/prof sync) ----
	graph_name = "pair"

	def graph_ready(self):
		'''True if initialized enough to be captured (kernels compiled, params up).'''
		return (self.cp is not None and self.kernel_force_write is not None
		        and self.kernel_force_add is not None and not self.params_changed)

	def graph_unready_reason(self):
		if self.cp is None or self.kernel_force_write is None or self.kernel_force_add is None:
			return "pair_not_initialized"
		if self.params_changed:
			return "pair_params_changed"
		return "unsupported"

	def graph_signature(self):
		'''Identity of device pointers/shape a captured graph bakes in; if it
		changes (neighbor-list rebuild/expand) the graph must be re-captured.'''
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		return (id(self.info.d_pos), id(self.info.d_force), id(self.d_params),
		        id(d_nlist), id(d_nsize), int(d_nlist.shape[1]),
		        int(self.info.ntypes), int(self.info.npa))

	def _graph_common_force_args(self):
		d_nlist = self.nlist.data.d_neighbor_list
		d_nsize = self.nlist.data.d_neighbor_size
		nlist_width = int(d_nlist.shape[1])
		pos   = self._dev("pos", self.info.d_pos)
		box   = self._dev("box", self.info.d_box)
		force = self._dev("force", self.info.d_force)
		nlist = d_nlist.ravel()       # CuPy-native nlist consumed directly
		nsize = d_nsize
		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		common = (
			np.int32(self.info.npa), pos, self.d_params, box,
			np.int32(self.info.ntypes), np.int32(self.nparams),
			nsize, nlist, np.int32(nlist_width), force,
		)
		return grid, block, common

	def graph_force_write(self):
		'''Launch ONLY the force-only overwrite pair kernel on the current cupy
		stream.  No nlist.calculate, no prof, no virial_potential.'''
		grid, block, common = self._graph_common_force_args()
		self.kernel_force_write(grid, block, common)

	def graph_force_add(self):
		'''Launch ONLY the force-only additive pair kernel on the current cupy
		stream.  Used by CUDA Graph multi-force graphs after captured zero_force_vp.'''
		grid, block, common = self._graph_common_force_args()
		self.kernel_force_add(grid, block, common)
