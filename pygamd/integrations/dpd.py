'''
CuPy-native minimal DPDGWVV integrator.

This is a lightweight Groot-Warren velocity-Verlet style DPD integrator for the
CuPy-native branch.  It follows the useful GALAMOST data-flow without copying
GALAMOST's class architecture:

    firststep:
        velo = vel_old + lambda_gw * dt * force_old / mass   # predictor velocity
        vel  = vel_old + 0.5       * dt * force_old / mass   # half-step velocity
        pos  = pos + vel * dt

    RawDPD force:
        reads CuInfo.d_velo when velocity_source='velo'

    secondstep:
        vel = vel + 0.5 * dt * force_new / mass

The d_velo array is a stable CuInfo core field, so RawDPD can use the same
RawKernel ABI under CUDA Graph capture/replay.  This file intentionally does not
add rigid body support or a parallel DPD integrator family.
'''

import numpy as np

from pygamd.profiler import prof
from pygamd import engine_metrics as _sp
from pygamd import engine_timing as _t
from pygamd import launch_trace as _lt
from pygamd.backend.cupy import import_cupy, as_cupy_array
from pygamd.backend.cupy.runtime import (
    DPDGWVV_DEFAULT_LAMBDA, VELOCITY_SOURCE_VELO,
    DeviceViewCache, launch_config_1d, require_cuinfo,
)


_DPDGWVV_BLOCK = 256

_KERNEL_SRC = r'''
extern "C" __global__ void dpdgwvv_first_step(
        const int nme, const int* __restrict__ member,
        float* __restrict__ pos, float* __restrict__ vel,
        float* __restrict__ velo,
        const float* __restrict__ force, int* __restrict__ image,
        const float* __restrict__ box, const float dt, const float lambda_gw)
{
    int t = blockIdx.x * blockDim.x + threadIdx.x;
    if (t < nme)
    {
        int idx = member[t];
        float m = vel[idx * 4 + 3];
        float ax = force[idx * 3 + 0] / m;
        float ay = force[idx * 3 + 1] / m;
        float az = force[idx * 3 + 2] / m;

        float voldx = vel[idx * 4 + 0];
        float voldy = vel[idx * 4 + 1];
        float voldz = vel[idx * 4 + 2];

        // Predictor velocity used by DPD dissipative force.
        velo[idx * 4 + 0] = voldx + lambda_gw * ax * dt;
        velo[idx * 4 + 1] = voldy + lambda_gw * ay * dt;
        velo[idx * 4 + 2] = voldz + lambda_gw * az * dt;
        velo[idx * 4 + 3] = m;

        float vx = voldx + 0.5f * ax * dt;
        float vy = voldy + 0.5f * ay * dt;
        float vz = voldz + 0.5f * az * dt;

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

extern "C" __global__ void dpdgwvv_second_step(
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


class gwvv:
    '''CuPy-native minimal DPDGWVV integrator.

    The constructor accepts an optional dpd_force wrapper/backend.  When provided,
    it switches that RawDPD force to velocity_source='velo'.  This keeps the user
    API clear while avoiding a parallel DPD integrator family.
    '''

    graph_name = "gwvv"

    def __init__(self, info, ci, dpd_force=None, lambda_gw=DPDGWVV_DEFAULT_LAMBDA):
        self.info = info
        self.cu = require_cuinfo(info, "integration.dpdgwvv")
        if not self.cu.has_field("velo"):
            raise RuntimeError("integration.dpdgwvv requires CuInfo.d_velo")
        self.ci = ci
        self.nme = int(ci.ps.nme)
        self.block_size = _DPDGWVV_BLOCK
        self.grid, self.block = launch_config_1d(self.nme, self.block_size)
        self.lambda_gw = float(lambda_gw)
        self.cp = None
        self.module = None
        self.k_first = None
        self.k_second = None
        self._wrap_cache = DeviceViewCache(as_cupy_array)
        self.dpd_force = None
        if dpd_force is not None:
            self.attach_dpd_force(dpd_force)

    def attach_dpd_force(self, dpd_force):
        backend = getattr(dpd_force, "data", dpd_force)
        if not hasattr(backend, "setVelocitySource"):
            raise TypeError("integration.dpdgwvv dpd_force must be pygamd.force.dpd or RawDPD backend")
        backend.setVelocitySource(VELOCITY_SOURCE_VELO)
        self.dpd_force = backend

    def setLambda(self, lambda_gw):
        self.lambda_gw = float(lambda_gw)

    def register(self, timestep):
        pass

    def _ensure(self):
        if self.cp is not None:
            return
        cp = import_cupy("RawDPDGWVV")
        self.cp = cp
        self.module = cp.RawModule(code=_KERNEL_SRC)
        self.k_first = self.module.get_function("dpdgwvv_first_step")
        self.k_second = self.module.get_function("dpdgwvv_second_step")

    def _wrap(self, key, arr):
        return self._wrap_cache.get(key, arr, ravel=True)

    def _dev(self, name, fallback):
        if self.cu is not None:
            return self.cu.dev(name)
        return self._wrap(name, fallback)

    def firststep(self, timestep):
        self._ensure()
        if _sp.ENABLED:
            _sp.inc("integration_firststep_count")
            _sp.inc("integration_dpdgwvv_firststep_count")
        member = self._wrap("member", self.ci.ps.d_member)
        pos = self._dev("pos", self.info.d_pos)
        vel = self._dev("vel", self.info.d_vel)
        velo = self._dev("velo", self.info.cu.d_velo)
        force = self._dev("force", self.info.d_force)
        image = self._dev("image", self.info.d_image)
        box = self._dev("box", self.info.d_box)
        dt = np.float32(self.info.dt)
        lam = np.float32(self.lambda_gw)
        grid, block = self.grid, self.block
        prof.begin("rawdpdgwvv_firststep")
        if _lt.LAUNCH: _lt.launch("integrator_dpdgwvv_firststep")
        if _t.ENABLED: _t.gpu_begin("integrator_dpdgwvv_firststep")
        self.k_first(grid, block, (np.int32(self.nme), member, pos, vel, velo, force, image, box, dt, lam))
        if _t.ENABLED: _t.gpu_end("integrator_dpdgwvv_firststep")
        prof.end("rawdpdgwvv_firststep")

    def secondstep(self, timestep):
        self._ensure()
        if _sp.ENABLED:
            _sp.inc("integration_secondstep_count")
            _sp.inc("integration_dpdgwvv_secondstep_count")
        member = self._wrap("member", self.ci.ps.d_member)
        vel = self._dev("vel", self.info.d_vel)
        force = self._dev("force", self.info.d_force)
        dt = np.float32(self.info.dt)
        grid, block = self.grid, self.block
        prof.begin("rawdpdgwvv_secondstep")
        if _lt.LAUNCH: _lt.launch("integrator_dpdgwvv_secondstep")
        if _t.ENABLED: _t.gpu_begin("integrator_dpdgwvv_secondstep")
        self.k_second(grid, block, (np.int32(self.nme), member, vel, force, dt))
        if _t.ENABLED: _t.gpu_end("integrator_dpdgwvv_secondstep")
        prof.end("rawdpdgwvv_secondstep")

    # --- CUDA Graph helpers --------------------------------------------------
    def graph_ready(self):
        try:
            self._ensure()
        except Exception:
            return False
        return self.k_first is not None and self.k_second is not None and self.cu.has_field("velo")

    def graph_unready_reason(self):
        if not self.cu.has_field("velo"):
            return "gwvv_velo_missing"
        if self.cp is None or self.k_first is None or self.k_second is None:
            return "gwvv_not_initialized"
        return "unsupported"

    def graph_signature(self):
        return ("gwvv", int(self.nme), int(self.block_size), float(self.info.dt),
                float(self.lambda_gw), id(self._dev("pos", self.info.d_pos)),
                id(self._dev("vel", self.info.d_vel)), id(self._dev("velo", self.info.cu.d_velo)),
                id(self._dev("force", self.info.d_force)), id(self._dev("image", self.info.d_image)),
                id(self._dev("box", self.info.d_box)), id(self._wrap("member", self.ci.ps.d_member)))

    def graph_firststep(self):
        self._ensure()
        member = self._wrap("member", self.ci.ps.d_member)
        pos = self._dev("pos", self.info.d_pos)
        vel = self._dev("vel", self.info.d_vel)
        velo = self._dev("velo", self.info.cu.d_velo)
        force = self._dev("force", self.info.d_force)
        image = self._dev("image", self.info.d_image)
        box = self._dev("box", self.info.d_box)
        dt = np.float32(self.info.dt)
        lam = np.float32(self.lambda_gw)
        grid, block = self.grid, self.block
        self.k_first(grid, block, (np.int32(self.nme), member, pos, vel, velo, force, image, box, dt, lam))

    def graph_secondstep(self):
        self._ensure()
        member = self._wrap("member", self.ci.ps.d_member)
        vel = self._dev("vel", self.info.d_vel)
        force = self._dev("force", self.info.d_force)
        dt = np.float32(self.info.dt)
        grid, block = self.grid, self.block
        self.k_second(grid, block, (np.int32(self.nme), member, vel, force, dt))
