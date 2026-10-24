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
CuPy-native dihedral force.

Consumes the CuInfo dihedral topology (get_dihedral_topology) and accumulates
into d_force / d_virial_potential with cupy.RawKernel + atomicAdd.  One thread
per dihedral (index-based connectivity i-j-k-l).  No numba.

Built-ins:
    periodic: U = k*(1 - cos(phi - phi0)), param [k, phi0_deg], proper only.
    harmonic: legacy pygamd-numba harmonic, param [k, t0_deg], proper/improper.

Custom device potentials are created with pygamd.potentials.dihedral_raw(...).  They
use the same ABI as pygamd-numba: func(cos_abcd, sin_abcd, param, fp).  For
custom dihedrals the backend appends term id to param internally:
    0 = proper, 1 = improper.
'''

import math
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp
from pygamd import engine_timing as _t
from pygamd import launch_trace as _lt
from pygamd.forces import potentials as _potentials
from pygamd.backend.cupy import import_cupy
from pygamd.backend.cupy.topology import get_dihedral_topology


_DIHEDRAL_KERNEL_TEMPLATE = r'''
%(device_code)s

extern "C" __global__ void %(kernel_name)s(
		const int    ndih,
		const float* __restrict__ pos,
		const int*   __restrict__ dlist,
		const int*   __restrict__ dtype,
		const float* __restrict__ params,
		const int    nparams,
		const float bx, const float by, const float bz,
		float* force, float* virial_potential)
{
	int d = blockIdx.x * blockDim.x + threadIdx.x;
	if (d < ndih)
	{
		int ia = dlist[d * 4 + 0];
		int ib = dlist[d * 4 + 1];
		int ic = dlist[d * 4 + 2];
		int id = dlist[d * 4 + 3];
		const float* pm = params + dtype[d] * nparams;
		int term = %(term_expr)s;

		float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;
		float ab0 = pos[ia*4+0]-pos[ib*4+0], ab1 = pos[ia*4+1]-pos[ib*4+1], ab2 = pos[ia*4+2]-pos[ib*4+2];
		float cb0 = pos[ic*4+0]-pos[ib*4+0], cb1 = pos[ic*4+1]-pos[ib*4+1], cb2 = pos[ic*4+2]-pos[ib*4+2];
		float dc0 = pos[id*4+0]-pos[ic*4+0], dc1 = pos[id*4+1]-pos[ic*4+1], dc2 = pos[id*4+2]-pos[ic*4+2];
		if      (ab0 >=  bxh) ab0 -= bx; else if (ab0 < -bxh) ab0 += bx;
		if      (ab1 >=  byh) ab1 -= by; else if (ab1 < -byh) ab1 += by;
		if      (ab2 >=  bzh) ab2 -= bz; else if (ab2 < -bzh) ab2 += bz;
		if      (cb0 >=  bxh) cb0 -= bx; else if (cb0 < -bxh) cb0 += bx;
		if      (cb1 >=  byh) cb1 -= by; else if (cb1 < -byh) cb1 += by;
		if      (cb2 >=  bzh) cb2 -= bz; else if (cb2 < -bzh) cb2 += bz;
		if      (dc0 >=  bxh) dc0 -= bx; else if (dc0 < -bxh) dc0 += bx;
		if      (dc1 >=  byh) dc1 -= by; else if (dc1 < -byh) dc1 += by;
		if      (dc2 >=  bzh) dc2 -= bz; else if (dc2 < -bzh) dc2 += bz;

		float fp[2];
		fp[0] = 0.0f;
		fp[1] = 0.0f;
		float ffa0 = 0.0f, ffa1 = 0.0f, ffa2 = 0.0f;
		float ffb0 = 0.0f, ffb1 = 0.0f, ffb2 = 0.0f;
		float ffc0 = 0.0f, ffc1 = 0.0f, ffc2 = 0.0f;
		float ffd0 = 0.0f, ffd1 = 0.0f, ffd2 = 0.0f;

		if (term == 0)
		{
			float cbm0 = -cb0, cbm1 = -cb1, cbm2 = -cb2;
			float aa0 = ab1*cbm2 - ab2*cbm1, aa1 = ab2*cbm0 - ab0*cbm2, aa2 = ab0*cbm1 - ab1*cbm0;
			float bb0 = dc1*cbm2 - dc2*cbm1, bb1 = dc2*cbm0 - dc0*cbm2, bb2 = dc0*cbm1 - dc1*cbm0;
			float aa_rsq = aa0*aa0 + aa1*aa1 + aa2*aa2;
			float bb_rsq = bb0*bb0 + bb1*bb1 + bb2*bb2;
			float g_rsq = cbm0*cbm0 + cbm1*cbm1 + cbm2*cbm2;
			float g_r = sqrtf(g_rsq);
			float g_rinv = (g_r > 0.0f) ? 1.0f / g_r : 0.0f;
			float aa_rsqinv = (aa_rsq > 0.0f) ? 1.0f / aa_rsq : 0.0f;
			float bb_rsqinv = (bb_rsq > 0.0f) ? 1.0f / bb_rsq : 0.0f;
			float ab_rinv = sqrtf(aa_rsqinv * bb_rsqinv);

			float cos_abcd = (aa0*bb0 + aa1*bb1 + aa2*bb2) * ab_rinv;
			float sin_abcd = g_r * ab_rinv * (aa0*dc0 + aa1*dc1 + aa2*dc2);
			if (cos_abcd >  1.0f) cos_abcd =  1.0f;
			if (cos_abcd < -1.0f) cos_abcd = -1.0f;

			%(proper_func)s(cos_abcd, sin_abcd, pm, fp);
			float f = fp[0];

			float fg = ab0*cbm0 + ab1*cbm1 + ab2*cbm2;
			float hg = dc0*cbm0 + dc1*cbm1 + dc2*cbm2;
			float fga = fg * aa_rsqinv * g_rinv;
			float hgb = hg * bb_rsqinv * g_rinv;
			float gaa = -aa_rsqinv * g_r;
			float gbb =  bb_rsqinv * g_r;

			float df0 = gaa*aa0, df1 = gaa*aa1, df2 = gaa*aa2;
			float dg0 = fga*aa0 - hgb*bb0, dg1 = fga*aa1 - hgb*bb1, dg2 = fga*aa2 - hgb*bb2;
			float dth0 = gbb*bb0, dth1 = gbb*bb1, dth2 = gbb*bb2;

			float ss0 = f*dg0, ss1 = f*dg1, ss2 = f*dg2;
			ffa0 = f*df0; ffa1 = f*df1; ffa2 = f*df2;
			ffb0 = ss0 - ffa0; ffb1 = ss1 - ffa1; ffb2 = ss2 - ffa2;
			ffd0 = f*dth0; ffd1 = f*dth1; ffd2 = f*dth2;
			ffc0 = -ss0 - ffd0; ffc1 = -ss1 - ffd1; ffc2 = -ss2 - ffd2;
		}
		else if (term == 1)
		{
			float r1 = rsqrtf(ab0*ab0 + ab1*ab1 + ab2*ab2);
			float r2 = rsqrtf(cb0*cb0 + cb1*cb1 + cb2*cb2);
			float r3 = rsqrtf(dc0*dc0 + dc1*dc1 + dc2*dc2);
			float ss1 = r1 * r1;
			float ss2 = r2 * r2;
			float ss3 = r3 * r3;

			float c0 =  (ab0*dc0 + ab1*dc1 + ab2*dc2) * r1 * r3;
			float c1 =  (ab0*cb0 + ab1*cb1 + ab2*cb2) * r1 * r2;
			float c2 = -(dc0*cb0 + dc1*cb1 + dc2*cb2) * r3 * r2;

			float s1 = 1.0f - c1*c1;
			if (s1 < 0.001f) s1 = 0.001f;
			s1 = 1.0f / s1;
			float s2 = 1.0f - c2*c2;
			if (s2 < 0.001f) s2 = 0.001f;
			s2 = 1.0f / s2;

			float s12 = sqrtf(s1 * s2);
			float cos_abcd = (c1*c2 + c0) * s12;
			if (cos_abcd >  1.0f) cos_abcd =  1.0f;
			if (cos_abcd < -1.0f) cos_abcd = -1.0f;
			float sin_abcd = sqrtf(1.0f - cos_abcd*cos_abcd);
			if (sin_abcd < 0.001f) sin_abcd = 0.001f;

			%(improper_func)s(cos_abcd, sin_abcd, pm, fp);
			float f = -fp[0] / sin_abcd;
			cos_abcd = cos_abcd * f;
			s12 = s12 * f;
			float a11 = cos_abcd * ss1 * s1;
			float a22 = -ss2 * (2.0f*c0*s12 - cos_abcd*(s1+s2));
			float a33 = cos_abcd * ss3 * s2;
			float a12 = -r1*r2 * (c1*cos_abcd*s1 + c2*s12);
			float a13 = -r1*r3 * s12;
			float a23 =  r2*r3 * (c2*cos_abcd*s2 + c1*s12);

			float ss0 = a22*cb0 + a23*dc0 + a12*ab0;
			float ssy = a22*cb1 + a23*dc1 + a12*ab1;
			float ssz = a22*cb2 + a23*dc2 + a12*ab2;

			ffa0 = a12*cb0 + a13*dc0 + a11*ab0;
			ffa1 = a12*cb1 + a13*dc1 + a11*ab1;
			ffa2 = a12*cb2 + a13*dc2 + a11*ab2;
			ffb0 = -ss0 - ffa0;
			ffb1 = -ssy - ffa1;
			ffb2 = -ssz - ffa2;
			ffd0 = a23*cb0 + a33*dc0 + a13*ab0;
			ffd1 = a23*cb1 + a33*dc1 + a13*ab1;
			ffd2 = a23*cb2 + a33*dc2 + a13*ab2;
			ffc0 = ss0 - ffd0;
			ffc1 = ssy - ffd1;
			ffc2 = ssz - ffd2;
		}
		else
		{
			return;
		}

		atomicAdd(&force[ia*3+0], ffa0); atomicAdd(&force[ia*3+1], ffa1); atomicAdd(&force[ia*3+2], ffa2);
		atomicAdd(&force[ib*3+0], ffb0); atomicAdd(&force[ib*3+1], ffb1); atomicAdd(&force[ib*3+2], ffb2);
		atomicAdd(&force[ic*3+0], ffc0); atomicAdd(&force[ic*3+1], ffc1); atomicAdd(&force[ic*3+2], ffc2);
		atomicAdd(&force[id*3+0], ffd0); atomicAdd(&force[id*3+1], ffd1); atomicAdd(&force[id*3+2], ffd2);

		float vx = ab0*ffa0 + cb0*ffc0 + (dc0 + cb0)*ffd0;
		float vy = ab1*ffa1 + cb1*ffc1 + (dc1 + cb1)*ffd1;
		float vz = ab2*ffa2 + cb2*ffc2 + (dc2 + cb2)*ffd2;
		float vir = (1.0f / 12.0f) * (vx + vy + vz);
		float ph  = 0.25f * fp[1];
		atomicAdd(&virial_potential[ia*2+0], vir); atomicAdd(&virial_potential[ia*2+1], ph);
		atomicAdd(&virial_potential[ib*2+0], vir); atomicAdd(&virial_potential[ib*2+1], ph);
		atomicAdd(&virial_potential[ic*2+0], vir); atomicAdd(&virial_potential[ic*2+1], ph);
		atomicAdd(&virial_potential[id*2+0], vir); atomicAdd(&virial_potential[id*2+1], ph);
	}
}
'''

def _func_name(func):
	if func is None:
		return "periodic"
	if isinstance(func, str):
		return func
	return getattr(func, "name", getattr(func, "__name__", getattr(func, "kind", "periodic")))


def _as_potential(func):
	if getattr(func, "kind", None) == "dihedral_raw":
		return func
	name = _func_name(func)
	if name in ("periodic", "harmonic"):
		return _potentials.get_dihedral_builtin(name)
	raise RuntimeError(
		"CuPy dihedral force supports built-in 'periodic'/'harmonic' or "
		"a custom pygamd.potentials.dihedral_raw(...); got %r." % (name,))


def _term_id(term):
	if term == "proper":
		return 0
	if term == "improper":
		return 1
	raise RuntimeError("Error, dihedral term %r; candidates are 'proper' and 'improper'." % (term,))


class dihedral:
	'''CuPy-native dihedral force backend.'''

	def __init__(self, info, func="periodic"):
		self.info = info
		self.potential = _as_potential(func)
		self.cos_factor = -1.0
		self.func = self.potential.name
		self.user_nparams = int(self.potential.nparams)
		self.param_style = getattr(self.potential, "param_style", None)
		self.proper_only = bool(getattr(self.potential, "proper_only", False))
		if self.param_style in ("dihedral_periodic", "dihedral_harmonic"):
			self.nparams = int(getattr(self.potential, "device_nparams", self.user_nparams))
		else:
			self.nparams = self.user_nparams + 1
		self.fast_math = bool(getattr(self.potential, "fast_math", False))
		self.device_code = self.potential.code
		self.proper_func = self.potential.name
		self.improper_func = self.potential.name
		self.kernel_name = "pygamd_dihedral_" + self.potential.name
		self.term_from_params = True

		self.topo = get_dihedral_topology(info)
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
		raise RuntimeError("CuPy dihedral: unknown dihedral type %r (known types: %s)" % (name, tm))

	def setCosFactor(self, factor):
		self.cos_factor = float(factor)
		if self.param_style in ("dihedral_periodic", "dihedral_harmonic"):
			idx = 3 if self.param_style == "dihedral_periodic" else 4
			for pm in self.params:
				if len(pm) > idx:
					pm[idx] = self.cos_factor
		self.params_changed = True

	def setParams(self, dihedral_type, param, term="proper", check=None):
		term_id = _term_id(term)
		tid = self._typeid(dihedral_type)

		if self.proper_only and term_id != 0:
			raise NotImplementedError(
				"CuPy dihedral '%s' is proper-only. Use harmonic or "
				"pygamd.potentials.dihedral_raw(...) for improper terms." % self.func)

		if self.param_style == "dihedral_periodic":
			if len(param) != self.user_nparams:
				raise RuntimeError("CuPy dihedral periodic expects 2 parameters [k, phi0_deg]; got %d." % len(param))
			k = float(param[0])
			phi0_rad = math.pi * float(param[1]) / 180.0
			self.params[tid] = [k, math.cos(phi0_rad), math.sin(phi0_rad),
			                    float(self.cos_factor), float(term_id)]
		elif self.param_style == "dihedral_harmonic":
			if len(param) != self.user_nparams:
				raise RuntimeError("CuPy dihedral harmonic expects 2 parameters [k, t0_deg]; got %d." % len(param))
			k = float(param[0])
			t0_rad = math.pi * float(param[1]) / 180.0
			self.params[tid] = [
				k, t0_rad, math.cos(t0_rad), math.sin(t0_rad),
				float(self.cos_factor), float(term_id)
			]
		else:
			if len(param) != self.user_nparams:
				raise RuntimeError("CuPy dihedral potential '%s' expects %d parameters; got %d."
				                   % (self.potential.name, self.user_nparams, len(param)))
			self.params[tid] = [float(x) for x in param] + [float(term_id)]

		if check is not None:
			check[tid] = True
		self.params_changed = True

	def _compile(self):
		term_expr = "(int)(pm[nparams - 1] + 0.5f)" if self.term_from_params else "0"
		source = _DIHEDRAL_KERNEL_TEMPLATE % {
			"device_code": self.device_code,
			"kernel_name": self.kernel_name,
			"term_expr": term_expr,
			"proper_func": self.proper_func,
			"improper_func": self.improper_func,
		}
		options = ('--use_fast_math',) if self.fast_math else ()
		try:
			kernel = self.cp.RawKernel(source, self.kernel_name, options=options)
			kernel.compile()
			return kernel
		except Exception as e:
			raise RuntimeError(
				"Failed to compile CuPy dihedral potential '%s'.\n"
				"Expected device signature:\n"
				"    __device__ void %s(float cos_abcd, float sin_abcd, const float* param, float* fp)\n"
				"Original error: %s" % (self.func, self.proper_func, repr(e))
			) from e

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("CuPy dihedral force")
		if self.kernel is None:
			self.kernel = self._compile()
		return self.cp

	def calculate(self, timestep):
		if self.topo.count == 0:
			return
		cp = self._ensure()
		prof.begin("dihedral_total")
		prof.count("dihedral_total")
		if _sp.ENABLED:
			_sp.inc("dihedral_compute_count")
		if self.params_changed:
			self.d_params = cp.asarray(np.asarray(self.params, dtype=np.float32).ravel())
			self.params_changed = False
		nd = int(self.topo.count)
		dlist = self.topo.d_indices.ravel()
		dtype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((nd + self.block_size - 1) // self.block_size,)
		if _lt.LAUNCH:
			_lt.launch("dihedral_force")
		if _t.ENABLED:
			_t.gpu_begin("dihedral_force")
		self.kernel(grid, (self.block_size,),
		            (np.int32(nd), pos, dlist, dtype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
		if _t.ENABLED:
			_t.gpu_end("dihedral_force")
		prof.end("dihedral_total")

	graph_name = "dihedral"

	def graph_ready(self):
		return self.topo.count == 0 or (self.cp is not None and self.kernel is not None and not self.params_changed)

	def graph_unready_reason(self):
		if self.topo.count == 0:
			return "none"
		if self.cp is None or self.kernel is None:
			return "dihedral_not_initialized"
		if self.params_changed:
			return "dihedral_params_changed"
		return "unsupported"

	def graph_signature(self):
		return ("dihedral", self.func, bool(self.fast_math), int(self.topo.count),
		        id(self.topo.d_indices), id(self.topo.d_typeid), id(self.d_params),
		        id(self.info.cu.dev("pos")), id(self.info.cu.dev("force")),
		        id(self.info.cu.dev("virial_potential")), tuple(float(x) for x in self.info.box),
		        int(self.nparams))

	def graph_force_add(self):
		if self.topo.count == 0:
			return
		nd = int(self.topo.count)
		dlist = self.topo.d_indices.ravel()
		dtype = self.topo.d_typeid
		pos = self.info.cu.dev("pos")
		force = self.info.cu.dev("force")
		vp = self.info.cu.dev("virial_potential")
		box = self.info.box
		grid = ((nd + self.block_size - 1) // self.block_size,)
		self.kernel(grid, (self.block_size,),
		            (np.int32(nd), pos, dlist, dtype, self.d_params, np.int32(self.nparams),
		             np.float32(box[0]), np.float32(box[1]), np.float32(box[2]), force, vp))
