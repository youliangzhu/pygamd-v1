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
CuPy-native PPPM reciprocal-space electrostatics.

This force is the long-range partner for nonbonded_c(func="lj_coulomb_erfc" or
"coulomb_erfc").  The real-space force supplies erfc(kappa*r)/r; this module
adds the smooth reciprocal term with a particle-particle particle-mesh solve:

  1. cubic B-spline charge assignment to a periodic mesh (RawKernel);
  2. rfftn(rho) with cuFFT;
  3. multiply by the Ewald Green function with B-spline deconvolution;
  4. irfftn(phi_k) to mesh potential;
  5. analytic-differentiation force gather from the mesh potential (RawKernel).

The module is additive: it accumulates into d_force and d_virial_potential after
the application-level zero-force step.  CUDA Graph capture is intentionally not
advertised because cuFFT workspaces and temporary arrays make a stable graph
contract much stricter than the pair/bond kernels.
'''

import math
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp
from pygamd import engine_timing as _t
from pygamd import launch_trace as _lt
from pygamd.backend.cupy import import_cupy
from pygamd.backend.cupy.runtime import require_cuinfo


_INV_SQRT_PI = 0.5641895835477563


def _next_fast_len(n):
	n = int(max(1, n))
	m = n
	while True:
		x = m
		for p in (2, 3, 5, 7):
			while x % p == 0:
				x //= p
		if x == 1:
			return m
		m += 1


def _normalize_grid(grid, box, mesh_spacing):
	if grid is not None:
		if isinstance(grid, int):
			g = (grid, grid, grid)
		else:
			g = tuple(int(x) for x in grid)
		if len(g) != 3 or min(g) < 4:
			raise RuntimeError("PPPM grid must be an int or a 3-tuple with all dimensions >= 4.")
		return g
	spacing = float(mesh_spacing)
	if spacing <= 0.0:
		raise RuntimeError("PPPM mesh_spacing must be positive when grid is not specified.")
	return tuple(_next_fast_len(max(4, int(math.ceil(float(box[i]) / spacing)))) for i in range(3))


def _collect_exclusion_pairs_host(info, exclusion):
	kinds = set(exclusion or ())
	if not kinds:
		return np.zeros((0, 2), dtype=np.int32)
	cu = getattr(info, "cu", None)
	topology = getattr(cu, "topology", {}) if cu is not None else {}
	pairs = set()

	def add_pair(a, b):
		a = int(a); b = int(b)
		if a == b:
			return
		if a > b:
			a, b = b, a
		pairs.add((a, b))

	if "bond" in kinds and "bond" in topology:
		for a, b in topology["bond"].host_tags[:, :2]:
			add_pair(a, b)
	if "angle" in kinds and "angle" in topology:
		for a, _b, c in topology["angle"].host_tags[:, :3]:
			add_pair(a, c)
	if "dihedral" in kinds and "dihedral" in topology:
		for a, _b, _c, d in topology["dihedral"].host_tags[:, :4]:
			add_pair(a, d)
	if not pairs:
		return np.zeros((0, 2), dtype=np.int32)
	return np.asarray(sorted(pairs), dtype=np.int32)


_PPPM_SRC = r'''
__device__ __forceinline__ int pppm_wrap(int i, int n)
{
	if (i >= n) i -= n;
	if (i < 0)  i += n;
	return i;
}

__device__ __forceinline__ void pppm_cubic_weights(float u, float* w, float* dw)
{
	float u2 = u * u;
	float u3 = u2 * u;
	float omu = 1.0f - u;
	float omu2 = omu * omu;

	w[0] = (omu2 * omu) * 0.16666666666666666f;
	w[1] = (3.0f*u3 - 6.0f*u2 + 4.0f) * 0.16666666666666666f;
	w[2] = (-3.0f*u3 + 3.0f*u2 + 3.0f*u + 1.0f) * 0.16666666666666666f;
	w[3] = u3 * 0.16666666666666666f;

	dw[0] = -0.5f * omu2;
	dw[1] = 1.5f*u2 - 2.0f*u;
	dw[2] = -1.5f*u2 + u + 0.5f;
	dw[3] = 0.5f*u2;
}

__device__ __forceinline__ void pppm_scaled_coord(
		float x, float box, int n, float* s, int* base, float* u)
{
	float nf = (float)n;
	float sf = (x + 0.5f * box) * nf / box;
	sf -= floorf(sf / nf) * nf;
	int cell = (int)floorf(sf);
	*s = sf;
	*base = cell - 1;
	*u = sf - (float)cell;
}

extern "C" __global__ void pppm_assign_cubic(
		const int npa,
		const float* __restrict__ pos,
		const float* __restrict__ charge,
		const float* __restrict__ box,
		const int nx, const int ny, const int nz,
		float* __restrict__ rho)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= npa) return;

	float sx, sy, sz, ux, uy, uz;
	int bx0, by0, bz0;
	pppm_scaled_coord(pos[i*4+0], box[0], nx, &sx, &bx0, &ux);
	pppm_scaled_coord(pos[i*4+1], box[1], ny, &sy, &by0, &uy);
	pppm_scaled_coord(pos[i*4+2], box[2], nz, &sz, &bz0, &uz);

	float wx[4], wy[4], wz[4], dxw[4], dyw[4], dzw[4];
	pppm_cubic_weights(ux, wx, dxw);
	pppm_cubic_weights(uy, wy, dyw);
	pppm_cubic_weights(uz, wz, dzw);

	float qi = charge[i];
	for (int a = 0; a < 4; a++)
	{
		int ix = pppm_wrap(bx0 + a, nx);
		for (int b = 0; b < 4; b++)
		{
			int iy = pppm_wrap(by0 + b, ny);
			float wxy = wx[a] * wy[b] * qi;
			for (int c = 0; c < 4; c++)
			{
				int iz = pppm_wrap(bz0 + c, nz);
				int gid = (ix * ny + iy) * nz + iz;
				atomicAdd(&rho[gid], wxy * wz[c]);
			}
		}
	}
}

extern "C" __global__ void pppm_gather_cubic(
		const int npa,
		const float* __restrict__ pos,
		const float* __restrict__ charge,
		const float* __restrict__ box,
		const int nx, const int ny, const int nz,
		const float* __restrict__ phi,
		const float coulomb_eff,
		const float kappa,
		const int include_energy,
		float* __restrict__ force,
		float* __restrict__ virial_potential)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= npa) return;

	float sx, sy, sz, ux, uy, uz;
	int bx0, by0, bz0;
	pppm_scaled_coord(pos[i*4+0], box[0], nx, &sx, &bx0, &ux);
	pppm_scaled_coord(pos[i*4+1], box[1], ny, &sy, &by0, &uy);
	pppm_scaled_coord(pos[i*4+2], box[2], nz, &sz, &bz0, &uz);

	float wx[4], wy[4], wz[4], dxw[4], dyw[4], dzw[4];
	pppm_cubic_weights(ux, wx, dxw);
	pppm_cubic_weights(uy, wy, dyw);
	pppm_cubic_weights(uz, wz, dzw);

	float inv_dx = ((float)nx) / box[0];
	float inv_dy = ((float)ny) / box[1];
	float inv_dz = ((float)nz) / box[2];
	float phi_i = 0.0f;
	float dphix = 0.0f, dphiy = 0.0f, dphiz = 0.0f;

	for (int a = 0; a < 4; a++)
	{
		int ix = pppm_wrap(bx0 + a, nx);
		for (int b = 0; b < 4; b++)
		{
			int iy = pppm_wrap(by0 + b, ny);
			for (int c = 0; c < 4; c++)
			{
				int iz = pppm_wrap(bz0 + c, nz);
				int gid = (ix * ny + iy) * nz + iz;
				float p = phi[gid];
				float wabc = wx[a] * wy[b] * wz[c];
				phi_i += wabc * p;
				dphix += dxw[a] * inv_dx * wy[b] * wz[c] * p;
				dphiy += wx[a] * dyw[b] * inv_dy * wz[c] * p;
				dphiz += wx[a] * wy[b] * dzw[c] * inv_dz * p;
			}
		}
	}

	float qi = charge[i];
	force[i*3+0] += -qi * dphix;
	force[i*3+1] += -qi * dphiy;
	force[i*3+2] += -qi * dphiz;

	if (include_energy)
	{
		const float inv_sqrt_pi = 0.5641895835477563f;
		float self_u = coulomb_eff * kappa * inv_sqrt_pi * qi * qi;
		virial_potential[i*2+1] += 0.5f * qi * phi_i - self_u;
	}
}

extern "C" __global__ void pppm_exclusion_long_correction(
		const int npairs,
		const int* __restrict__ tag_pairs,
		const int* __restrict__ rtag,
		const float* __restrict__ pos,
		const float* __restrict__ charge,
		const float* __restrict__ box,
		const float coulomb_eff,
		const float kappa,
		const int include_energy,
		float* __restrict__ force,
		float* __restrict__ virial_potential)
{
	int p = blockIdx.x * blockDim.x + threadIdx.x;
	if (p >= npairs) return;

	int tag_i = tag_pairs[p*2+0];
	int tag_j = tag_pairs[p*2+1];
	int i = rtag[tag_i];
	int j = rtag[tag_j];

	float dx = pos[i*4+0] - pos[j*4+0];
	float dy = pos[i*4+1] - pos[j*4+1];
	float dz = pos[i*4+2] - pos[j*4+2];
	float bx = box[0], by = box[1], bz = box[2];
	float bxh = bx * 0.5f, byh = by * 0.5f, bzh = bz * 0.5f;
	if      (dx >=  bxh) dx -= bx; else if (dx < -bxh) dx += bx;
	if      (dy >=  byh) dy -= by; else if (dy < -byh) dy += by;
	if      (dz >=  bzh) dz -= bz; else if (dz < -bzh) dz += bz;

	float rsq = dx*dx + dy*dy + dz*dz;
	if (rsq <= 0.0f) return;
	float r2inv = 1.0f / rsq;
	float rinv = sqrtf(r2inv);
	float r = rsq * rinv;
	float qiqj = charge[i] * charge[j];
	float kr = kappa * r;
	const float inv_sqrt_pi = 0.5641895835477563f;
	float erf_kr = erff(kr);
	float exp_kr2 = expf(-kr * kr);

	// Subtract the smooth long-range pair that the mesh includes for this
	// excluded bonded pair: -coulomb_eff*q_i*q_j*erf(kappa*r)/r.
	float fcoef = -coulomb_eff * qiqj *
		(erf_kr * r2inv * rinv - 2.0f * kappa * inv_sqrt_pi * exp_kr2 * r2inv);
	float fx = fcoef * dx, fy = fcoef * dy, fz = fcoef * dz;
	atomicAdd(&force[i*3+0],  fx);
	atomicAdd(&force[i*3+1],  fy);
	atomicAdd(&force[i*3+2],  fz);
	atomicAdd(&force[j*3+0], -fx);
	atomicAdd(&force[j*3+1], -fy);
	atomicAdd(&force[j*3+2], -fz);

	if (include_energy)
	{
		float u = -coulomb_eff * qiqj * erf_kr * rinv;
		float vir = (1.0f / 6.0f) * rsq * fcoef;
		float ph = 0.5f * u;
		atomicAdd(&virial_potential[i*2+0], vir);
		atomicAdd(&virial_potential[j*2+0], vir);
		atomicAdd(&virial_potential[i*2+1], ph);
		atomicAdd(&virial_potential[j*2+1], ph);
	}
}
'''


class pppm:
	'''PPPM reciprocal electrostatic force backend.'''

	graph_name = "pppm"

	def __init__(self, info, grid=None, mesh_spacing=1.0, kappa=None, rcut=None,
	             epsilonr=1.0, coulomb_const=138.935, exclusion=None,
	             block_size=128, deconvolve=True):
		self.info = info
		self.cu = require_cuinfo(info, "PPPM electrostatics")
		self.grid_request = grid
		self.mesh_spacing = float(mesh_spacing)
		if kappa is None:
			if rcut is None:
				raise RuntimeError("PPPM requires kappa or rcut.  If kappa is omitted, it uses 3.0/rcut.")
			kappa = 3.0 / float(rcut)
		self.kappa = float(kappa)
		if self.kappa <= 0.0:
			raise RuntimeError("PPPM kappa must be positive.")
		self.epsilonr = float(epsilonr)
		if self.epsilonr == 0.0:
			raise RuntimeError("PPPM epsilonr must be non-zero.")
		self.coulomb_eff = float(coulomb_const) / self.epsilonr
		self.exclusion = tuple(exclusion) if exclusion else None
		self.block_size = int(block_size)
		if self.block_size < 32 or self.block_size > 1024 or self.block_size % 32 != 0:
			self.block_size = 128
		self.deconvolve = bool(deconvolve)
		self.can_fuse_zero_force = False

		self.cp = None
		self.mod = None
		self.k_assign = None
		self.k_gather = None
		self.k_excl = None
		self.grid = None
		self._last_box = None
		self.d_rho = None
		self.d_phi = None
		self.d_green = None
		self.d_excl_pairs = None
		self.excl_count = 0
		self.nblocks = math.ceil(int(self.info.npa) / self.block_size)
		self._excl_pairs_host = _collect_exclusion_pairs_host(info, self.exclusion)

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("PPPM electrostatics")
			self.mod = self.cp.RawModule(code=_PPPM_SRC, options=('--use_fast_math',))
			self.k_assign = self.mod.get_function("pppm_assign_cubic")
			self.k_gather = self.mod.get_function("pppm_gather_cubic")
			self.k_excl = self.mod.get_function("pppm_exclusion_long_correction")
			if self._excl_pairs_host.shape[0] > 0:
				self.d_excl_pairs = self.cp.asarray(self._excl_pairs_host, dtype=self.cp.int32)
				self.excl_count = int(self._excl_pairs_host.shape[0])
		return self.cp

	def _sinc(self, x):
		cp = self.cp
		ax = cp.abs(x)
		return cp.where(ax < 1.0e-7, cp.ones_like(x), cp.sin(x) / x)

	def _build_green(self, box, grid):
		cp = self.cp
		nx, ny, nz = grid
		lx, ly, lz = box
		vol = float(lx * ly * lz)
		mesh_count = float(nx * ny * nz)
		kx = (2.0 * math.pi * cp.fft.fftfreq(nx, d=lx / nx)).astype(cp.float32)
		ky = (2.0 * math.pi * cp.fft.fftfreq(ny, d=ly / ny)).astype(cp.float32)
		kz = (2.0 * math.pi * cp.fft.rfftfreq(nz, d=lz / nz)).astype(cp.float32)
		k2 = (kx[:, None, None] * kx[:, None, None]
		      + ky[None, :, None] * ky[None, :, None]
		      + kz[None, None, :] * kz[None, None, :])
		mask = k2 > 0.0
		green = cp.zeros(k2.shape, dtype=cp.float32)
		base = (mesh_count / vol) * (4.0 * math.pi * self.coulomb_eff)
		k2_safe = cp.where(mask, k2, cp.ones_like(k2))
		val = base * cp.exp(-k2_safe / np.float32(4.0 * self.kappa * self.kappa)) / k2_safe
		if self.deconvolve:
			wx = self._sinc(0.5 * kx * np.float32(lx / nx))[:, None, None]
			wy = self._sinc(0.5 * ky * np.float32(ly / ny))[None, :, None]
			wz = self._sinc(0.5 * kz * np.float32(lz / nz))[None, None, :]
			window2 = (wx * wy * wz) ** np.float32(8.0)
			val = val / cp.maximum(window2, np.float32(1.0e-8))
		green[mask] = val[mask]
		green[0, 0, 0] = np.float32(0.0)
		return green

	def _ensure_mesh(self):
		cp = self._ensure()
		box = tuple(float(x) for x in self.info.box[:3])
		grid = _normalize_grid(self.grid_request, box, self.mesh_spacing)
		if self.d_rho is not None and self.grid == grid and self._last_box == box:
			return
		self.grid = grid
		self._last_box = box
		self.d_rho = cp.zeros(grid, dtype=cp.float32)
		self.d_phi = cp.zeros(grid, dtype=cp.float32)
		self.d_green = self._build_green(box, grid)
		if _sp.ENABLED:
			_sp.set_meta("pppm_grid", "%d,%d,%d" % grid)
			_sp.set_meta("pppm_kappa", float(self.kappa))
			_sp.set_meta("pppm_mesh_spacing", float(self.mesh_spacing))
			_sp.set_meta("pppm_exclusion_pairs", int(self.excl_count))

	def _need_energy(self):
		cp = getattr(self.info, "compute_properties", None)
		if not isinstance(cp, dict):
			return True
		return bool(
			cp.get("potential", False)
			or cp.get("pressure", False)
			or cp.get("virial", False)
			or cp.get("stress_tensor", False)
		)

	def set_fuse_zero_force(self, flag):
		pass

	def calculate(self, timestep):
		if getattr(self.info, "d_charge", None) is None:
			raise RuntimeError("PPPM electrostatics requires per-particle charges in the snapshot.")
		cp = self._ensure()
		self._ensure_mesh()
		prof.begin("pppm_total")
		prof.count("pppm_total")
		if _sp.ENABLED:
			_sp.inc("pppm_compute_count")
		nx, ny, nz = self.grid
		grid = (int(self.nblocks),)
		block = (int(self.block_size),)
		pos = self.cu.dev("pos")
		charge = self.cu.dev("charge")
		box = self.cu.dev("box")
		force = self.cu.dev("force")
		vp = self.cu.dev("virial_potential")

		self.d_rho.fill(0.0)
		if _lt.LAUNCH: _lt.launch("pppm_assign")
		if _t.ENABLED: _t.gpu_begin("pppm_assign")
		self.k_assign(grid, block,
		              (np.int32(self.info.npa), pos, charge, box,
		               np.int32(nx), np.int32(ny), np.int32(nz), self.d_rho.ravel()))
		if _t.ENABLED: _t.gpu_end("pppm_assign")

		if _lt.LAUNCH: _lt.launch("pppm_fft")
		if _t.ENABLED: _t.gpu_begin("pppm_fft")
		rho_k = cp.fft.rfftn(self.d_rho)
		phi_k = rho_k * self.d_green
		self.d_phi[...] = cp.fft.irfftn(phi_k, s=self.grid)
		if _t.ENABLED: _t.gpu_end("pppm_fft")

		include_energy = np.int32(1 if self._need_energy() else 0)
		if _lt.LAUNCH: _lt.launch("pppm_gather")
		if _t.ENABLED: _t.gpu_begin("pppm_gather")
		self.k_gather(grid, block,
		              (np.int32(self.info.npa), pos, charge, box,
		               np.int32(nx), np.int32(ny), np.int32(nz), self.d_phi.ravel(),
		               np.float32(self.coulomb_eff), np.float32(self.kappa),
		               include_energy, force, vp))
		if _t.ENABLED: _t.gpu_end("pppm_gather")

		if self.excl_count > 0:
			eg = ((self.excl_count + self.block_size - 1) // self.block_size,)
			if _lt.LAUNCH: _lt.launch("pppm_exclusion_correction")
			if _t.ENABLED: _t.gpu_begin("pppm_exclusion_correction")
			self.k_excl(eg, block,
			            (np.int32(self.excl_count), self.d_excl_pairs.ravel(),
			             self.cu.dev("rtag"), pos, charge, box,
			             np.float32(self.coulomb_eff), np.float32(self.kappa),
			             include_energy, force, vp))
			if _t.ENABLED: _t.gpu_end("pppm_exclusion_correction")
		prof.end("pppm_total")

	def graph_ready(self):
		return False

	def graph_unready_reason(self):
		return "pppm_fft_eager_only"

	def graph_signature(self):
		return ("pppm", int(self.info.npa), self.grid, float(self.kappa), float(self.epsilonr))
