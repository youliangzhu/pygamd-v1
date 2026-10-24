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
CuPy-native angle force.

Consumes the CuInfo angle topology (get_angle_topology) and accumulates into
d_force / d_virial_potential with a cupy.RawKernel + atomicAdd.  One thread per
angle (index-based connectivity i-j-k, j = center).  No numba.

Built-in 'harmonic':  U(theta) = 0.5 * k * (theta - theta0)^2
    param = [k, theta0_degrees]   (theta0 is given in DEGREES, converted to
    radians internally -- matches legacy pygamd).  Force geometry and the
    sin(theta) ~ 0 guard (minimum 0.001) follow the legacy numba angle kernel.

Custom device potentials are created with pygamd.potentials.angle_raw(...).  They use
the same ABI as pygamd-numba: func(cos_abc, sin_abc, param, fp).

virial_potential layout (npa, 2): [:,0]=virial, [:,1]=potential.  Per angle, each
of the three particles receives (1/6)*(d_ab.fab + d_cb.fcb) of virial and (1/3)*U
of potential (totals (1/2)*(...) and U) -- identical to the legacy per-particle
angle accumulation.

count == 0 -> safe no-op.  cupy imported lazily; importing needs no CUDA driver.
'''

import math
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.forces import potentials as _potentials
from pygamd.backend.cupy import import_cupy
from pygamd.backend.cupy.topology import get_angle_topology


_ANGLE_KERNEL_TEMPLATE = r'''
%(device_code)s

extern "C" __global__ void %(kernel_name)s(
		const int    nangles,
		const float* __restrict__ pos,
		const int*   __restrict__ alist,
		const int*   __restrict__ atype,
		const float* __restrict__ params,
		const int    nparams,
		const float bx, const float by, const float bz,
		float* force, float* virial_potential)
{
	int a = blockIdx.x * blockDim.x + threadIdx.x;
	if (a < nangles)
	{
		int i = alist[a * 3 + 0];     // a
		int j = alist[a * 3 + 1];     // b (center)
		int kk = alist[a * 3 + 2];    // c
		const float* pm = params + atype[a] * nparams;

		float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;
		// d_ab = pa - pb (i - j),  d_cb = pc - pb (k - j)
		float ab0 = pos[i*4+0]-pos[j*4+0], ab1 = pos[i*4+1]-pos[j*4+1], ab2 = pos[i*4+2]-pos[j*4+2];
		float cb0 = pos[kk*4+0]-pos[j*4+0], cb1 = pos[kk*4+1]-pos[j*4+1], cb2 = pos[kk*4+2]-pos[j*4+2];
		if      (ab0 >=  bxh) ab0 -= bx; else if (ab0 < -bxh) ab0 += bx;
		if      (ab1 >=  byh) ab1 -= by; else if (ab1 < -byh) ab1 += by;
		if      (ab2 >=  bzh) ab2 -= bz; else if (ab2 < -bzh) ab2 += bz;
		if      (cb0 >=  bxh) cb0 -= bx; else if (cb0 < -bxh) cb0 += bx;
		if      (cb1 >=  byh) cb1 -= by; else if (cb1 < -byh) cb1 += by;
		if      (cb2 >=  bzh) cb2 -= bz; else if (cb2 < -bzh) cb2 += bz;

		float rsq_ab = ab0*ab0 + ab1*ab1 + ab2*ab2;
		float rsq_cb = cb0*cb0 + cb1*cb1 + cb2*cb2;
		float r_ab = sqrtf(rsq_ab), r_cb = sqrtf(rsq_cb);

		float c = (ab0*cb0 + ab1*cb1 + ab2*cb2) / (r_ab * r_cb);
		if (c >  1.0f) c =  1.0f;
		if (c < -1.0f) c = -1.0f;
		float s = sqrtf(1.0f - c * c);
		if (s < 0.001f) s = 0.001f;
		float sinv = 1.0f / s;

		float fp[2];
		fp[0] = 0.0f;
		fp[1] = 0.0f;
		%(func_name)s(c, s, pm, fp);

		float aa  = -fp[0] * sinv;
		float a11 = aa * c / rsq_ab;
		float a12 = -aa / (r_ab * r_cb);
		float a22 = aa * c / rsq_cb;

		float fab0 = a11*ab0 + a12*cb0, fab1 = a11*ab1 + a12*cb1, fab2 = a11*ab2 + a12*cb2;
		float fcb0 = a22*cb0 + a12*ab0, fcb1 = a22*cb1 + a12*ab1, fcb2 = a22*cb2 + a12*ab2;

		// forces: i(a) += fab ; k(c) += fcb ; j(b) -= (fab+fcb)
		atomicAdd(&force[i*3+0], fab0);  atomicAdd(&force[i*3+1], fab1);  atomicAdd(&force[i*3+2], fab2);
		atomicAdd(&force[kk*3+0], fcb0); atomicAdd(&force[kk*3+1], fcb1); atomicAdd(&force[kk*3+2], fcb2);
		atomicAdd(&force[j*3+0], -(fab0+fcb0)); atomicAdd(&force[j*3+1], -(fab1+fcb1)); atomicAdd(&force[j*3+2], -(fab2+fcb2));

		float V = (ab0*fab0 + cb0*fcb0) + (ab1*fab1 + cb1*fcb1) + (ab2*fab2 + cb2*fcb2);
		float vir = (1.0f / 6.0f) * V;   // per-particle share (x3 -> V/2 total)
		float ph  = (1.0f / 3.0f) * fp[1]; // per-particle share (x3 -> U total)
		atomicAdd(&virial_potential[i*2+0], vir);  atomicAdd(&virial_potential[i*2+1], ph);
		atomicAdd(&virial_potential[j*2+0], vir);  atomicAdd(&virial_potential[j*2+1], ph);
		atomicAdd(&virial_potential[kk*2+0], vir); atomicAdd(&virial_potential[kk*2+1], ph);
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
	if getattr(func, "kind", None) == "angle_raw":
		return func
	name = _func_name(func)
	if name in ("harmonic", "harmonic_cos"):
		return _potentials.get_angle_builtin(name)
	raise RuntimeError(
		"CuPy angle force supports built-in 'harmonic'/'harmonic_cos' or a "
		"custom pygamd.potentials.angle_raw(...); got %r." % (name,))


class angle:
	'''CuPy-native angle force backend.'''

	def __init__(self, info, func="harmonic"):
		self.info = info
		self.potential = _as_potential(func)
		self.func = self.potential.name
		self.nparams = int(self.potential.nparams)
		self.fast_math = bool(getattr(self.potential, "fast_math", False))
		self.param_style = getattr(self.potential, "param_style", None)
		self.device_code = self.potential.code
		self.func_name = self.potential.name
		self.kernel_name = "pygamd_angle_" + self.potential.name
		self.topo = get_angle_topology(info)
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
		raise RuntimeError("RawAngle: unknown angle type %r (known types: %s)" % (name, tm))

	def setParams(self, angle_type, param, check=None):
		tid = self._typeid(angle_type)
		if len(param) != self.nparams:
			if self.param_style in ("angle_harmonic", "angle_harmonic_cos"):
				raise RuntimeError("CuPy angle %s expects 2 parameters [k, theta0_deg]; got %d."
				                   % (self.func, len(param)))
			raise RuntimeError("CuPy angle potential '%s' expects %d parameters; got %d."
			                   % (self.potential.name, self.nparams, len(param)))
		if self.param_style in ("angle_harmonic", "angle_harmonic_cos"):
			k = float(param[0])
			t0_rad = math.pi * float(param[1]) / 180.0
			if self.param_style == "angle_harmonic_cos":
				self.params[tid] = [k, math.cos(t0_rad)]
			else:
				self.params[tid] = [k, t0_rad]
		else:
			self.params[tid] = [float(x) for x in param]
		if check is not None:
			check[tid] = True
		self.params_changed = True

	def _compile(self):
		source = _ANGLE_KERNEL_TEMPLATE % {
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
				"Failed to compile CuPy angle potential '%s'.\n"
				"Expected device signature:\n"
				"    __device__ void %s(float cos_abc, float sin_abc, const float* param, float* fp)\n"
				"Original error: %s" % (self.func, self.func_name, repr(e))
			) from e

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("RawAngle force")
		if self.kernel is None:
			self.kernel = self._compile()
		return self.cp

	def calculate(self, timestep):
		if self.topo.count == 0:
			return
		cp = self._ensure()
		prof.begin("angle_total")
		prof.count("angle_total")
		if _sp.ENABLED: _sp.inc("angle_compute_count")
		if self.params_changed:
			self.d_params = cp.asarray(np.asarray(self.params, dtype=np.float32).ravel())
			self.params_changed = False
		na = int(self.topo.count)
		alist = self.topo.d_indices.ravel()
		atype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((na + self.block_size - 1) // self.block_size,)
		if _lt.LAUNCH: _lt.launch("angle_force")
		if _t.ENABLED: _t.gpu_begin("angle_force")
		self.kernel(grid, (self.block_size,),
		            (np.int32(na), pos, alist, atype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
		if _t.ENABLED: _t.gpu_end("angle_force")
		prof.end("angle_total")

	# --- CUDA Graph CUDA-graph helpers (additive bonded force segment) ----------
	graph_name = "angle"

	def graph_ready(self):
		return self.topo.count == 0 or (self.cp is not None and self.kernel is not None and not self.params_changed)

	def graph_unready_reason(self):
		if self.topo.count == 0:
			return "none"
		if self.cp is None or self.kernel is None:
			return "angle_not_initialized"
		if self.params_changed:
			return "angle_params_changed"
		return "unsupported"

	def graph_signature(self):
		return ("angle", self.func, bool(self.fast_math), int(self.topo.count), id(self.topo.d_indices), id(self.topo.d_typeid),
		        id(self.d_params), id(self.info.cu.dev("pos")), id(self.info.cu.dev("force")),
		        id(self.info.cu.dev("virial_potential")), tuple(float(x) for x in self.info.box),
		        int(self.nparams))

	def graph_force_add(self):
		if self.topo.count == 0:
			return
		na = int(self.topo.count)
		alist = self.topo.d_indices.ravel()
		atype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((na + self.block_size - 1) // self.block_size,)
		self.kernel(grid, (self.block_size,),
		            (np.int32(na), pos, alist, atype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
