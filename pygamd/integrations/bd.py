'''
CuPy-native BD (Brownian / Langevin) integrator.

Velocity-Verlet position update + a Langevin secondstep that adds a per-particle
dissipative (-gamma*v) and random force (Saru RNG, deterministic from
seed+timestep+tag) before the half-kick.  Mirrors the legacy numba BD:
    coeff = sqrt(6 * gamma * T / dt)
    rf = rand[-1,1] * coeff - gamma * v
    v += 0.5 * (force + rf) / m * dt
Per-type friction `gamma` via setParams(type, gamma).  Random numbers use the
Saru RNG keyed by (seed+timestep, tag, offset) -- tag-based so it is stable under
a future Sort.  cupy.RawKernel; no numba.

Variant temperature: setT(value) accepts a float or a VariantLinear-like object.
The graph path captures stable device scalar pointers for both seed and temperature;
graph_prepare_replay(timestep) updates them before graph.launch().
'''

import math
import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.backend.cupy import import_cupy, as_cupy_array
from pygamd.backend.cupy.rng import SARU_SRC


_BD_BLOCK = 256

_KERNEL_SRC = SARU_SRC + r'''
extern "C" __global__ void bd_first_step(
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
		float ax = force[idx * 3 + 0] / m, ay = force[idx * 3 + 1] / m, az = force[idx * 3 + 2] / m;
		float vx = vel[idx * 4 + 0], vy = vel[idx * 4 + 1], vz = vel[idx * 4 + 2];
		float px = pos[idx * 4 + 0] + vx * dt + 0.5f * ax * dt * dt;
		float py = pos[idx * 4 + 1] + vy * dt + 0.5f * ay * dt * dt;
		float pz = pos[idx * 4 + 2] + vz * dt + 0.5f * az * dt * dt;
		vx += 0.5f * ax * dt; vy += 0.5f * ay * dt; vz += 0.5f * az * dt;
		float bx = box[0], by = box[1], bz = box[2];
		int ix = image[idx * 3 + 0], iy = image[idx * 3 + 1], iz = image[idx * 3 + 2];
		if      (px >=  bx * 0.5f) { px -= bx; ix += 1; }
		else if (px <  -bx * 0.5f) { px += bx; ix -= 1; }
		if      (py >=  by * 0.5f) { py -= by; iy += 1; }
		else if (py <  -by * 0.5f) { py += by; iy -= 1; }
		if      (pz >=  bz * 0.5f) { pz -= bz; iz += 1; }
		else if (pz <  -bz * 0.5f) { pz += bz; iz -= 1; }
		pos[idx*4+0]=px; pos[idx*4+1]=py; pos[idx*4+2]=pz;
		vel[idx*4+0]=vx; vel[idx*4+1]=vy; vel[idx*4+2]=vz;
		image[idx*3+0]=ix; image[idx*3+1]=iy; image[idx*3+2]=iz;
	}
}

extern "C" __global__ void bd_second_step(
		const int nme, const int* __restrict__ member,
		const float* __restrict__ pos, float* __restrict__ vel,
		const float* __restrict__ force, const int* __restrict__ tag,
		const float* __restrict__ params,
		const unsigned int seed, const float T, const int dimension, const float dt)
{
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t < nme)
	{
		int idx = member[t];
		float m = vel[idx * 4 + 3];
		int   ti = (int)pos[idx * 4 + 3];
		float gamma = params[ti];
		float coeff = sqrtf(6.0f * gamma * T / dt);
		unsigned int tg = (unsigned int)tag[idx];

		Saru sx, sy, sz;
		saru3(&sx, seed,     tg,       seed + 1u);
		saru3(&sy, seed,     seed + 2u, tg);
		saru3(&sz, seed + 3u, tg,       seed);
		float rfx = saru_f(&sx, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 0];
		float rfy = saru_f(&sy, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 1];
		float rfz = saru_f(&sz, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 2];
		if (dimension == 2) rfz = 0.0f;

		float fx = force[idx * 3 + 0] + rfx;
		float fy = force[idx * 3 + 1] + rfy;
		float fz = force[idx * 3 + 2] + rfz;
		float minv = 1.0f / m;
		vel[idx * 4 + 0] += 0.5f * fx * minv * dt;
		vel[idx * 4 + 1] += 0.5f * fy * minv * dt;
		vel[idx * 4 + 2] += 0.5f * fz * minv * dt;
	}
}

extern "C" __global__ void bd_second_step_seedptr(
		const int nme, const int* __restrict__ member,
		const float* __restrict__ pos, float* __restrict__ vel,
		const float* __restrict__ force, const int* __restrict__ tag,
		const float* __restrict__ params,
		const unsigned int* __restrict__ seed_ptr,
		const float* __restrict__ T_ptr, const int dimension, const float dt)
{
	float T = T_ptr[0];
	int t = blockIdx.x * blockDim.x + threadIdx.x;
	if (t < nme)
	{
		int idx = member[t];
		float m = vel[idx * 4 + 3];
		int   ti = (int)pos[idx * 4 + 3];
		float gamma = params[ti];
		float coeff = sqrtf(6.0f * gamma * T / dt);
		unsigned int seed = seed_ptr[0];
		unsigned int tg = (unsigned int)tag[idx];

		Saru sx, sy, sz;
		saru3(&sx, seed,     tg,       seed + 1u);
		saru3(&sy, seed,     seed + 2u, tg);
		saru3(&sz, seed + 3u, tg,       seed);
		float rfx = saru_f(&sx, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 0];
		float rfy = saru_f(&sy, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 1];
		float rfz = saru_f(&sz, -1.0f, 1.0f) * coeff - gamma * vel[idx * 4 + 2];
		if (dimension == 2) rfz = 0.0f;

		float fx = force[idx * 3 + 0] + rfx;
		float fy = force[idx * 3 + 1] + rfy;
		float fz = force[idx * 3 + 2] + rfz;
		float minv = 1.0f / m;
		vel[idx * 4 + 0] += 0.5f * fx * minv * dt;
		vel[idx * 4 + 1] += 0.5f * fy * minv * dt;
		vel[idx * 4 + 2] += 0.5f * fz * minv * dt;
	}
}

'''


class bd:
	'''CuPy-native Brownian dynamics (Langevin).'''

	def __init__(self, info, ci, temp, seed=12345):
		self.info = info
		self.cu = getattr(info, "cu", None)
		self.ci = ci
		self.temp = temp
		self.seed = int(seed)
		self.nme = int(ci.ps.nme)
		self.block_size = _BD_BLOCK
		self.nblocks = math.ceil(self.nme / self.block_size)
		self.ntypes = int(info.ntypes)
		self.params = [1.0 for _ in range(self.ntypes)]   # per-type gamma (default 1)
		self.params_changed = True
		self.cp = None
		self.k_first = None
		self.k_second = None
		self.k_second_graph = None
		self.d_params = None
		self.d_graph_seed = None
		self.d_graph_temp = None
		self._wrap_cache = {}

	def setParams(self, typ, gamma):
		tid = self.info.convert_name_to_id(typ)
		self.params[tid] = float(gamma)
		self.params_changed = True

	def setT(self, temperature):
		self.temp = temperature

	def _get_T(self, timestep):
		if hasattr(self.temp, "getValue"):
			return float(self.temp.getValue(timestep))
		return float(self.temp)

	def register(self, timestep):
		pass

	def _ensure(self):
		if self.cp is None:
			self.cp = import_cupy("RawBD")
			mod = self.cp.RawModule(code=_KERNEL_SRC)
			self.k_first = mod.get_function("bd_first_step")
			self.k_second = mod.get_function("bd_second_step")
			self.k_second_graph = mod.get_function("bd_second_step_seedptr")

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
		prof.begin("rawbd_firststep")
		if _lt.LAUNCH: _lt.launch("integrator_bd")
		if _t.ENABLED: _t.gpu_begin("integrator_bd")
		self.k_first(grid, block, (np.int32(self.nme), member, pos, vel, force, image, box, dt))
		if _t.ENABLED: _t.gpu_end("integrator_bd")
		prof.end("rawbd_firststep")

	def secondstep(self, timestep):
		self._ensure()
		if _sp.ENABLED: _sp.inc("integration_secondstep_count")
		if self.params_changed:
			self.d_params = self.cp.asarray(np.asarray(self.params, dtype=np.float32))
			self.params_changed = False
		member = self._wrap("member", self.ci.ps.d_member)
		pos = self._dev("pos", self.info.d_pos)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		tag = self._dev("tag", self.info.d_tag)
		dt = np.float32(self.info.dt)
		T = np.float32(self._get_T(timestep))
		seed = np.uint32((self.seed + int(timestep)) & 0xffffffff)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		prof.begin("rawbd_secondstep")
		if _lt.LAUNCH: _lt.launch("integrator_bd")
		if _t.ENABLED: _t.gpu_begin("integrator_bd")
		self.k_second(grid, block, (np.int32(self.nme), member, pos, vel, force, tag,
		                            self.d_params, seed, T, np.int32(self.info.dimension), dt))
		if _t.ENABLED: _t.gpu_end("integrator_bd")
		prof.end("rawbd_secondstep")


	# --- engine metrics graph helpers -------------------------------------------------
	# BD uses a stochastic secondstep.  The eager path passes seed+timestep by value;
	# graph replay would otherwise freeze that scalar at capture time.  The graph
	# path therefore captures a device seed pointer and updates it before every
	# replay via graph_prepare_replay(), matching the eager RNG key.
	graph_name = "bd"

	def _temp_is_variant(self):
		return hasattr(self.temp, "getValue")

	def _ensure_params_device(self):
		self._ensure()
		if self.params_changed or self.d_params is None:
			self.d_params = self.cp.asarray(np.asarray(self.params, dtype=np.float32))
			self.params_changed = False
		return self.d_params

	def _ensure_graph_seed(self):
		self._ensure()
		if self.d_graph_seed is None:
			self.d_graph_seed = self.cp.zeros((1,), dtype=self.cp.uint32)
		return self.d_graph_seed

	def _ensure_graph_temp(self):
		self._ensure()
		if self.d_graph_temp is None:
			self.d_graph_temp = self.cp.zeros((1,), dtype=self.cp.float32)
		return self.d_graph_temp

	def graph_ready(self):
		try:
			self._ensure()
			self._ensure_params_device()
		except Exception:
			return False
		self._ensure_graph_seed()
		self._ensure_graph_temp()
		return (self.k_first is not None and self.k_second_graph is not None
		        and self.d_params is not None and self.d_graph_seed is not None
		        and self.d_graph_temp is not None)

	def graph_unready_reason(self):
		if self.cp is None or self.k_first is None or self.k_second_graph is None:
			return "bd_not_initialized"
		if self.params_changed or self.d_params is None:
			return "bd_params_changed"
		return "unsupported"

	def graph_target_value(self):
		# Stable signature target only.  The graph reads d_graph_temp[0], updated
		# by graph_prepare_replay(timestep), so VariantLinear does not force fallback.
		return float(self._get_T(0))

	def graph_prepare_replay(self, timestep):
		d_seed = self._ensure_graph_seed()
		d_temp = self._ensure_graph_temp()
		seed = np.uint32((int(self.seed) + int(timestep)) & 0xffffffff)
		T = np.float32(self._get_T(timestep))
		d_seed.fill(seed)
		d_temp.fill(T)
		if _sp.ENABLED:
			_sp.inc("cuda_graph_bd_dynamic_seed_update_count")
			_sp.inc("cuda_graph_temperature_scalar_update_count")
			_sp.inc("cuda_graph_bd_temperature_scalar_update_count")
			if self._temp_is_variant():
				_sp.inc("cuda_graph_variant_temperature_update_count")
				_sp.inc("cuda_graph_bd_variant_temperature_update_count")

	def graph_signature(self):
		temp_mode = "variant_temperature" if self._temp_is_variant() else "static_temperature"
		temp_sig = id(self.temp) if self._temp_is_variant() else float(self.graph_target_value())
		return ("bd", int(self.nme), int(self.block_size), float(self.info.dt),
		        temp_mode, temp_sig, float(self.graph_target_value()), int(self.seed), int(self.info.dimension),
		        id(self._dev("pos", self.info.d_pos)), id(self._dev("vel", self.info.d_vel)),
		        id(self._dev("force", self.info.d_force)), id(self._dev("image", self.info.d_image)),
		        id(self._dev("box", self.info.d_box)), id(self._dev("tag", self.info.d_tag)),
		        id(self._wrap("member", self.ci.ps.d_member)), id(self._ensure_params_device()),
		        id(self._ensure_graph_seed()), id(self._ensure_graph_temp()))

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
		self._ensure_params_device()
		self._ensure_graph_seed()
		self._ensure_graph_temp()
		member = self._wrap("member", self.ci.ps.d_member)
		pos = self._dev("pos", self.info.d_pos)
		vel = self._dev("vel", self.info.d_vel)
		force = self._dev("force", self.info.d_force)
		tag = self._dev("tag", self.info.d_tag)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks),); block = (int(self.block_size),)
		self.k_second_graph(grid, block, (np.int32(self.nme), member, pos, vel, force, tag,
		                                 self.d_params, self.d_graph_seed, self.d_graph_temp,
		                                 np.int32(self.info.dimension), dt))
