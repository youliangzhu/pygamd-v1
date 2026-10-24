'''
CuPy-native bond force.

Consumes the CuInfo bond topology (pygamd.backend.cupy.topology.get_bond_topology)
and accumulates into d_force / d_virial_potential with a cupy.RawKernel + atomicAdd.
One thread per bond (index-based connectivity).  No numba.

Built-in 'harmonic':  U(r) = 0.5 * k * (r - r0)^2 ,  param = [k, r0]
    fp0 = k * (r0/r - 1)         (force/r factor; F_on_i = fp0 * (r_i - r_j))
    matches the legacy numba bond ABI exactly (fp[0]=force/r, fp[1]=potential).

Custom device potentials are created with pygamd.potentials.bond_raw(...).  They use
the same ABI as pygamd-numba: func(rsq, param, fp).

virial_potential layout is the engine's (npa, 2): [:,0]=virial, [:,1]=potential.
Per bond, each of the two particles receives (1/6)*rsq*fp0 of virial and 0.5*U of
potential (totals (1/3)*rsq*fp0 and U) -- identical to the legacy per-particle
bond accumulation, so compute_properties / dump.data are unchanged.

count == 0  ->  safe no-op.  cupy is imported lazily; importing this module needs
no CUDA driver.
'''

import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.forces import potentials as _potentials
from pygamd.backend.cupy import import_cupy
from pygamd.backend.cupy.topology import get_bond_topology


_BOND_KERNEL_TEMPLATE = r'''
%(device_code)s

extern "C" __global__ void %(kernel_name)s(
		const int    nbonds,
		const float* __restrict__ pos,
		const int*   __restrict__ blist,
		const int*   __restrict__ btype,
		const float* __restrict__ params,
		const int    nparams,
		const float bx, const float by, const float bz,
		float* force, float* virial_potential)
{
	int b = blockIdx.x * blockDim.x + threadIdx.x;
	if (b < nbonds)
	{
		int i = blist[b * 2 + 0];
		int j = blist[b * 2 + 1];
		const float* pm = params + btype[b] * nparams;

		float dx = pos[i * 4 + 0] - pos[j * 4 + 0];
		float dy = pos[i * 4 + 1] - pos[j * 4 + 1];
		float dz = pos[i * 4 + 2] - pos[j * 4 + 2];
		float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;
		if      (dx >=  bxh) dx -= bx; else if (dx < -bxh) dx += bx;
		if      (dy >=  byh) dy -= by; else if (dy < -byh) dy += by;
		if      (dz >=  bzh) dz -= bz; else if (dz < -bzh) dz += bz;

		float rsq = dx * dx + dy * dy + dz * dz;
		float r = sqrtf(rsq);
		if (r < 1e-6f) return;

		float fp[2];
		fp[0] = 0.0f;
		fp[1] = 0.0f;
		%(func_name)s(rsq, pm, fp);

		float fx = fp[0] * dx, fy = fp[0] * dy, fz = fp[0] * dz;
		atomicAdd(&force[i * 3 + 0],  fx); atomicAdd(&force[i * 3 + 1],  fy); atomicAdd(&force[i * 3 + 2],  fz);
		atomicAdd(&force[j * 3 + 0], -fx); atomicAdd(&force[j * 3 + 1], -fy); atomicAdd(&force[j * 3 + 2], -fz);

		float vir = (1.0f / 6.0f) * rsq * fp[0]; // per-particle share (x2 -> 1/3 total)
		float ph  = 0.5f * fp[1];                // per-particle share (x2 -> U total)
		atomicAdd(&virial_potential[i * 2 + 0], vir); atomicAdd(&virial_potential[i * 2 + 1], ph);
		atomicAdd(&virial_potential[j * 2 + 0], vir); atomicAdd(&virial_potential[j * 2 + 1], ph);
	}
}
'''


def _func_name(func):
	if func is None:
		return "harmonic"
	if isinstance(func, str):
		return func
	return getattr(func, "name", getattr(func, "__name__", getattr(func, "kind", "harmonic")))


def _as_potential(func):
	if getattr(func, "kind", None) == "bond_raw":
		return func
	name = _func_name(func)
	if name == "harmonic":
		return _potentials.get_bond_builtin(name)
	raise RuntimeError(
		"CuPy bond force supports built-in 'harmonic' or a custom "
		"pygamd.potentials.bond_raw(...); got %r." % (name,))


class bond:
	'''CuPy-native bond force backend.'''

	def __init__(self, info, func="harmonic"):
		self.info = info
		self.potential = _as_potential(func)
		self.func = self.potential.name
		self.nparams = int(self.potential.nparams)
		self.fast_math = bool(getattr(self.potential, "fast_math", False))
		self.param_style = getattr(self.potential, "param_style", None)
		self.device_code = self.potential.code
		self.func_name = self.potential.name
		self.kernel_name = "pygamd_bond_" + self.potential.name
		self.topo = get_bond_topology(info)
		self.ntypes = max(len(self.topo.typemap), 1)
		self.params = [[0.0] * self.nparams for _ in range(self.ntypes)]
		self.params_changed = True
		self.block_size = 128
		self.cp = None
		self.kernel = None
		self.d_params = None

	def _typeid(self, name):
		tm = self.topo.typemap
		if name in tm:
			return tm.index(name)
		raise RuntimeError("RawBond: unknown bond type %r (known types: %s)" % (name, tm))

	def setParams(self, bond_type, param, check=None):
		tid = self._typeid(bond_type)
		if len(param) != self.nparams:
			if self.param_style == "bond_harmonic":
				raise RuntimeError("CuPy bond harmonic expects 2 parameters [k, r0]; got %d." % len(param))
			raise RuntimeError("CuPy bond potential '%s' expects %d parameters; got %d."
			                   % (self.potential.name, self.nparams, len(param)))
		self.params[tid] = [float(x) for x in param]
		if check is not None:
			check[tid] = True
		self.params_changed = True

	def _compile(self):
		source = _BOND_KERNEL_TEMPLATE % {
			"device_code": self.device_code,
			"kernel_name": self.kernel_name,
			"func_name": self.func_name,
		}
		options = ('--use_fast_math',) if self.fast_math else ()
		try:
			kernel = self.cp.RawKernel(source, self.kernel_name, options=options)
			kernel.compile()
			return kernel
		except Exception as e:
			raise RuntimeError(
				"Failed to compile CuPy bond potential '%s'.\n"
				"Expected device signature:\n"
				"    __device__ void %s(float rsq, const float* param, float* fp)\n"
				"Original error: %s" % (self.func, self.func_name, repr(e))
			) from e

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("RawBond force")
		if self.kernel is None:
			self.kernel = self._compile()
		return self.cp

	def calculate(self, timestep):
		if self.topo.count == 0:
			return  # safe no-op when there are no bonds
		cp = self._ensure()
		prof.begin("bond_total")
		prof.count("bond_total")
		if _sp.ENABLED: _sp.inc("bond_compute_count")
		if self.params_changed:
			self.d_params = cp.asarray(np.asarray(self.params, dtype=np.float32).ravel())
			self.params_changed = False
		nb = int(self.topo.count)
		blist = self.topo.d_indices.ravel()
		btype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((nb + self.block_size - 1) // self.block_size,)
		if _lt.LAUNCH: _lt.launch("bond_force")
		if _t.ENABLED: _t.gpu_begin("bond_force")
		self.kernel(grid, (self.block_size,),
		            (np.int32(nb), pos, blist, btype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
		if _t.ENABLED: _t.gpu_end("bond_force")
		prof.end("bond_total")

	# --- CUDA Graph CUDA-graph helpers (additive bonded force segment) ----------
	graph_name = "bond"

	def graph_ready(self):
		return self.topo.count == 0 or (self.cp is not None and self.kernel is not None and not self.params_changed)

	def graph_unready_reason(self):
		if self.topo.count == 0:
			return "none"
		if self.cp is None or self.kernel is None:
			return "bond_not_initialized"
		if self.params_changed:
			return "bond_params_changed"
		return "unsupported"

	def graph_signature(self):
		return ("bond", self.func, bool(self.fast_math), int(self.topo.count), id(self.topo.d_indices), id(self.topo.d_typeid),
		        id(self.d_params), id(self.info.cu.dev("pos")), id(self.info.cu.dev("force")),
		        id(self.info.cu.dev("virial_potential")), tuple(float(x) for x in self.info.box),
		        int(self.nparams))

	def graph_force_add(self):
		if self.topo.count == 0:
			return
		nb = int(self.topo.count)
		blist = self.topo.d_indices.ravel()
		btype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((nb + self.block_size - 1) // self.block_size,)
		self.kernel(grid, (self.block_size,),
		            (np.int32(nb), pos, blist, btype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
