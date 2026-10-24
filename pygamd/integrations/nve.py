'''
CuPy-native RawNVE integrator (velocity-Verlet, no thermostat).

Same lifecycle as RawNVT (firststep / [force compute by application] / secondstep)
but with no Nose-Hoover variable: it is RawNVT with xi = 0.  Consumes the CuInfo
CuPy arrays (d_pos/d_vel/d_force/d_image/d_box; mass = vel[:,3]) for the group's
members.  cupy.RawKernel; no numba.

    firststep : a = force/m ; v += 0.5*a*dt ; x += v*dt ; PBC + image update
    secondstep: a = force/m ; v += 0.5*a*dt
'''

import math
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.backend.cupy import import_cupy, as_cupy_array


_NVE_BLOCK = 256

_KERNEL_SRC = r'''
extern "C" __global__ void nve_first_step(
		const int nme, const int* __restrict__ member,
		float* __restrict__ pos, float* __restrict__ vel,
		const float* __restrict__ force, int* __restrict__ image,
		const float* __restrict__ box, const float dt)
{
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t < nme)
	{
		int idx = member[t];
		float m = vel[idx * 4 + 3];
		float ax = force[idx * 3 + 0] / m;
		float ay = force[idx * 3 + 1] / m;
		float az = force[idx * 3 + 2] / m;
		float vx = vel[idx * 4 + 0] + 0.5f * ax * dt;
		float vy = vel[idx * 4 + 1] + 0.5f * ay * dt;
		float vz = vel[idx * 4 + 2] + 0.5f * az * dt;
		float px = pos[idx * 4 + 0] + vx * dt;
		float py = pos[idx * 4 + 1] + vy * dt;
		float pz = pos[idx * 4 + 2] + vz * dt;
		float bx = box[0], by = box[1], bz = box[2];
		int ix = image[idx * 3 + 0], iy = image[idx * 3 + 1], iz = image[idx * 3 + 2];
		if      (px >=  bx * 0.5f) { px -= bx; ix += 1; }
		else if (px <  -bx * 0.5f) { px += bx; ix -= 1; }
		if      (py >=  by * 0.5f) { py -= by; iy += 1; }
		else if (py <  -by * 0.5f) { py += by; iy -= 1; }
		if      (pz >=  bz * 0.5f) { pz -= bz; iz += 1; }
		else if (pz <  -bz * 0.5f) { pz += bz; iz -= 1; }
		pos[idx * 4 + 0] = px; pos[idx * 4 + 1] = py; pos[idx * 4 + 2] = pz;
		vel[idx * 4 + 0] = vx; vel[idx * 4 + 1] = vy; vel[idx * 4 + 2] = vz;
		image[idx * 3 + 0] = ix; image[idx * 3 + 1] = iy; image[idx * 3 + 2] = iz;
	}
}

extern "C" __global__ void nve_second_step(
		const int nme, const int* __restrict__ member,
		float* __restrict__ vel, const float* __restrict__ force, const float dt)
{
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t < nme)
	{
		int idx = member[t];
		float m = vel[idx * 4 + 3];
		vel[idx * 4 + 0] += 0.5f * (force[idx * 3 + 0] / m) * dt;
		vel[idx * 4 + 1] += 0.5f * (force[idx * 3 + 1] / m) * dt;
		vel[idx * 4 + 2] += 0.5f * (force[idx * 3 + 2] / m) * dt;
	}
}
'''


class nve:
	'''CuPy-native velocity-Verlet NVE (same ctor shape as the legacy nve(info, ci)).'''

	def __init__(self, info, ci):
		self.info = info
		self.cu = getattr(info, "cu", None)
		self.ci = ci
		self.nme = int(ci.ps.nme)
		self.block_size = _NVE_BLOCK
		self.nblocks = math.ceil(self.nme / self.block_size)
		self.cp = None
		self.k_first = None
		self.k_second = None
		self._wrap_cache = {}

	def setT(self, temperature):
		# NVE has no thermostat; accept-and-ignore would be silent, so be explicit.
		raise RuntimeError("integration.nve has no temperature; use integration.nvt for a thermostat.")

	def register(self, timestep):
		pass

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("RawNVE")
			mod = self.cp.RawModule(code=_KERNEL_SRC)
			self.k_first = mod.get_function("nve_first_step")
			self.k_second = mod.get_function("nve_second_step")

	def _wrap(self, key, arr):
		c = self._wrap_cache.get(key)
		if c is not None and c[0] is arr:
			return c[1]
		w = as_cupy_array(arr, ravel=True)
		self._wrap_cache[key] = (arr, w)
		return w

	def _dev(self, name, fallback):
		if self.cu is not None:
			return self.cu.dev(name)
		return self._wrap(name, fallback)

	def firststep(self, timestep):
		self._ensure()
		if _sp.ENABLED: _sp.inc("integration_firststep_count")
		member = self._wrap("member", self.ci.ps.d_member)
		pos = self._dev("pos", self.info.d_pos)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		image = self._dev("image", self.info.d_image)
		box = self._dev("box", self.info.d_box)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		prof.begin("rawnve_firststep")
		if _lt.LAUNCH: _lt.launch("integrator_nve")
		if _t.ENABLED: _t.gpu_begin("integrator_nve")
		self.k_first(grid, block, (np.int32(self.nme), member, pos, vel, force, image, box, dt))
		if _t.ENABLED: _t.gpu_end("integrator_nve")
		prof.end("rawnve_firststep")

	def secondstep(self, timestep):
		self._ensure()
		if _sp.ENABLED: _sp.inc("integration_secondstep_count")
		member = self._wrap("member", self.ci.ps.d_member)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		prof.begin("rawnve_secondstep")
		if _lt.LAUNCH: _lt.launch("integrator_nve")
		if _t.ENABLED: _t.gpu_begin("integrator_nve")
		self.k_second(grid, block, (np.int32(self.nme), member, vel, force, dt))
		if _t.ENABLED: _t.gpu_end("integrator_nve")
		prof.end("rawnve_secondstep")


	# --- engine metrics graph helpers -------------------------------------------------
	# NVE has no thermostat stage.  The CUDA-graph engine therefore captures
	# firststep + force segment + secondstep directly.  These helpers launch the
	# same RawKernels as the eager path, with no profiler/timing side effects.
	graph_name = "nve"

	def graph_ready(self):
		try:
			self._ensure()
		except Exception:
			return False
		return self.k_first is not None and self.k_second is not None

	def graph_unready_reason(self):
		if self.cp is None or self.k_first is None or self.k_second is None:
			return "nve_not_initialized"
		return "unsupported"

	def graph_signature(self):
		return ("nve", int(self.nme), int(self.block_size), float(self.info.dt),
		        id(self._dev("pos", self.info.d_pos)), id(self._dev("vel", self.info.d_vel)),
		        id(self._dev("force", self.info.d_force)), id(self._dev("image", self.info.d_image)),
		        id(self._dev("box", self.info.d_box)), id(self._wrap("member", self.ci.ps.d_member)))

	def graph_firststep(self):
		self._ensure()
		member = self._wrap("member", self.ci.ps.d_member)
		pos = self._dev("pos", self.info.d_pos)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		image = self._dev("image", self.info.d_image)
		box = self._dev("box", self.info.d_box)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		self.k_first(grid, block, (np.int32(self.nme), member, pos, vel, force, image, box, dt))

	def graph_secondstep(self):
		self._ensure()
		member = self._wrap("member", self.ci.ps.d_member)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		self.k_second(grid, block, (np.int32(self.nme), member, vel, force, dt))
