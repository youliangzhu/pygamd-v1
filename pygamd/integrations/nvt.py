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
# RawKernel-native Nose-Hoover NVT integrator.
#
# Replaces the per-step host round-trip of the Numba NVT path:
#   * firststep / secondstep are cupy.RawKernels,
#   * the group temperature is reduced on the device (two-stage reduction),
#   * the thermostat variable xi lives in a device scalar d_xi[0] and is updated
#     by a device kernel; no device->host copy per step, no ci.calculate() per step.
#   * force is overwritten with force / mass inside the velocity steps, matching
#     the legacy numba NVT side effect that subsequent firststep calls observe.
#
# Physics is identical to integrations/nvt.py (Nose-Hoover):
#   firststep : a = force/m;  dn_inv = 1/(1 + dt/2*xi)
#               v = (v + dt/2*a) * dn_inv ;  x += v*dt ; PBC + image update
#   secondstep: temp = sum_members m*(vx^2+vy^2+vz^2) / nfreedom
#               xi += dt/tau^2 * (temp/target_temp - 1)
#               v = v + dt/2*(a - xi*v)
#
# cupy is imported lazily; `import pygamd` never depends on cupy.  Graph replay
# handles static and variant target temperatures through a stable one-element
# device scalar updated by graph_prepare_replay(), matching the DPD dynamic-seed
# replay mechanism instead of baking a host scalar into the captured graph.
# -----------------------------------------------------------------------------

import os
import math
import numpy as np
from pygamd.profiler import prof
from pygamd import engine_metrics as _sp   # engine metrics opt-in counters
from pygamd import engine_timing as _t     # engine timing opt-in timing decomposition
from pygamd import launch_trace as _lt   # launch trace opt-in launch trace
from pygamd.backend.cupy import import_cupy, as_cupy_array, asnumpy


_NVT_BLOCK = 256          # threads per block for the member kernels / reduction
_FINALIZE_THREADS = 256   # single-block reduction width for xi update
_FUSE_THREADS = 1024      # single-block width for the fused xi-update + secondstep


_KERNEL_SRC = r'''
extern "C" __global__ void nvt_first_step(
		const int npa_unused, const int nme, const int* __restrict__ member,
		float* __restrict__ pos, float* __restrict__ vel,
		float* __restrict__ force, int* __restrict__ image,
		const float* __restrict__ box, const float* __restrict__ xi, const float dt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nme)
	{
		int idx = member[i];
		float m  = vel[idx * 4 + 3];
		float ax = force[idx * 3 + 0] / m;
		float ay = force[idx * 3 + 1] / m;
		float az = force[idx * 3 + 2] / m;
		force[idx * 3 + 0] = ax;
		force[idx * 3 + 1] = ay;
		force[idx * 3 + 2] = az;

		float dn_inv = 1.0f / (1.0f + 0.5f * dt * xi[0]);
		float vx = (vel[idx * 4 + 0] + 0.5f * ax * dt) * dn_inv;
		float vy = (vel[idx * 4 + 1] + 0.5f * ay * dt) * dn_inv;
		float vz = (vel[idx * 4 + 2] + 0.5f * az * dt) * dn_inv;

		float px = pos[idx * 4 + 0] + vx * dt;
		float py = pos[idx * 4 + 1] + vy * dt;
		float pz = pos[idx * 4 + 2] + vz * dt;

		float bx = box[0], by = box[1], bz = box[2];
		int ix = image[idx * 3 + 0];
		int iy = image[idx * 3 + 1];
		int iz = image[idx * 3 + 2];

		if      (px >=  bx * 0.5f) { px -= bx; ix += 1; }
		else if (px <  -bx * 0.5f) { px += bx; ix -= 1; }
		if      (py >=  by * 0.5f) { py -= by; iy += 1; }
		else if (py <  -by * 0.5f) { py += by; iy -= 1; }
		if      (pz >=  bz * 0.5f) { pz -= bz; iz += 1; }
		else if (pz <  -bz * 0.5f) { pz += bz; iz -= 1; }

		pos[idx * 4 + 0] = px;  pos[idx * 4 + 1] = py;  pos[idx * 4 + 2] = pz;
		vel[idx * 4 + 0] = vx;  vel[idx * 4 + 1] = vy;  vel[idx * 4 + 2] = vz;
		image[idx * 3 + 0] = ix; image[idx * 3 + 1] = iy; image[idx * 3 + 2] = iz;
	}
}

extern "C" __global__ void nvt_ke_partial(
		const int nme, const int* __restrict__ member,
		const float* __restrict__ vel, float* __restrict__ partial)
{
	__shared__ float sm[256];
	int tx = threadIdx.x;
	int i  = blockIdx.x * blockDim.x + tx;
	float t = 0.0f;
	if (i < nme)
	{
		int idx = member[i];
		float vx = vel[idx * 4 + 0];
		float vy = vel[idx * 4 + 1];
		float vz = vel[idx * 4 + 2];
		float m  = vel[idx * 4 + 3];
		t = m * (vx * vx + vy * vy + vz * vz);
	}
	sm[tx] = t;
	__syncthreads();
	for (int offs = blockDim.x >> 1; offs > 0; offs >>= 1) {
		if (tx < offs) sm[tx] += sm[tx + offs];
		__syncthreads();
	}
	if (tx == 0) partial[blockIdx.x] = sm[0];
}

extern "C" __global__ void nvt_finalize_xi(
		const float* __restrict__ partial, const int nblocks, float* __restrict__ xi,
		const float dt_over_tau2, const float inv_nfreedom, const float target_temp)
{
	__shared__ float sm[256];
	int tx = threadIdx.x;
	float s = 0.0f;
	for (int o = tx; o < nblocks; o += blockDim.x) s += partial[o];
	sm[tx] = s;
	__syncthreads();
	for (int offs = blockDim.x >> 1; offs > 0; offs >>= 1) {
		if (tx < offs) sm[tx] += sm[tx + offs];
		__syncthreads();
	}
	if (tx == 0) {
		float temp = sm[0] * inv_nfreedom;
		xi[0] += dt_over_tau2 * (temp / target_temp - 1.0f);
		xi[1] = temp;   // debug cache: last device-side temperature
	}
}

extern "C" __global__ void nvt_finalize_xi_tptr(
		const float* __restrict__ partial, const int nblocks, float* __restrict__ xi,
		const float dt_over_tau2, const float inv_nfreedom,
		const float* __restrict__ target_temp_ptr)
{
	__shared__ float sm[256];
	int tx = threadIdx.x;
	float s = 0.0f;
	for (int o = tx; o < nblocks; o += blockDim.x) s += partial[o];
	sm[tx] = s;
	__syncthreads();
	for (int offs = blockDim.x >> 1; offs > 0; offs >>= 1) {
		if (tx < offs) sm[tx] += sm[tx + offs];
		__syncthreads();
	}
	if (tx == 0) {
		float target_temp = target_temp_ptr[0];
		float temp = sm[0] * inv_nfreedom;
		xi[0] += dt_over_tau2 * (temp / target_temp - 1.0f);
		xi[1] = temp;   // debug cache: last device-side temperature
	}
}

extern "C" __global__ void nvt_second_step(
		const int nme, const int* __restrict__ member,
		float* __restrict__ vel, float* __restrict__ force,
		const float* __restrict__ xi, const float dt)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i < nme)
	{
		int idx = member[i];
		float m  = vel[idx * 4 + 3];
		float ax = force[idx * 3 + 0] / m;
		float ay = force[idx * 3 + 1] / m;
		float az = force[idx * 3 + 2] / m;
		force[idx * 3 + 0] = ax;
		force[idx * 3 + 1] = ay;
		force[idx * 3 + 2] = az;
		float x  = xi[0];
		float vx = vel[idx * 4 + 0];
		float vy = vel[idx * 4 + 1];
		float vz = vel[idx * 4 + 2];
		vx = vx + 0.5f * (ax - x * vx) * dt;
		vy = vy + 0.5f * (ay - x * vy) * dt;
		vz = vz + 0.5f * (az - x * vz) * dt;
		vel[idx * 4 + 0] = vx;
		vel[idx * 4 + 1] = vy;
		vel[idx * 4 + 2] = vz;
	}
}

// Fused kernel: reduce KE partials -> update xi[0] -> (single-block
// barrier) -> apply secondstep to all members with the NEW xi.  Launched with a
// SINGLE block of _FUSE_THREADS threads so __syncthreads() validly orders the
// xi update before every member's velocity update (no cross-block sync).
// Physics is identical to nvt_finalize_xi followed by nvt_second_step.
extern "C" __global__ void nvt_finalize_xi_secondstep(
		const float* __restrict__ partial, const int nblocks,
		float* __restrict__ xi, const float dt_over_tau2,
		const float inv_nfreedom, const float target_temp,
		const int nme, const int* __restrict__ member,
		float* __restrict__ vel, float* __restrict__ force, const float dt)
{
	__shared__ float sm[1024];
	__shared__ float sh_xi;
	int tx = threadIdx.x;

	// stage 1: reduce KE partial sums (single block)
	float s = 0.0f;
	for (int o = tx; o < nblocks; o += blockDim.x) s += partial[o];
	sm[tx] = s;
	__syncthreads();
	for (int offs = blockDim.x >> 1; offs > 0; offs >>= 1) {
		if (tx < offs) sm[tx] += sm[tx + offs];
		__syncthreads();
	}
	// stage 2: thread 0 updates the thermostat variable
	if (tx == 0) {
		float temp = sm[0] * inv_nfreedom;
		float xnew = xi[0] + dt_over_tau2 * (temp / target_temp - 1.0f);
		xi[0] = xnew;
		xi[1] = temp;       // debug cache: last device-side temperature
		sh_xi = xnew;
	}
	__syncthreads();        // ensure new xi visible to all threads before use
	float x = sh_xi;

	// stage 3: secondstep velocity update for every member, using the NEW xi
	for (int i = tx; i < nme; i += blockDim.x) {
		int idx = member[i];
		float m  = vel[idx * 4 + 3];
		float ax = force[idx * 3 + 0] / m;
		float ay = force[idx * 3 + 1] / m;
		float az = force[idx * 3 + 2] / m;
		force[idx * 3 + 0] = ax;
		force[idx * 3 + 1] = ay;
		force[idx * 3 + 2] = az;
		float vx = vel[idx * 4 + 0];
		float vy = vel[idx * 4 + 1];
		float vz = vel[idx * 4 + 2];
		vx = vx + 0.5f * (ax - x * vx) * dt;
		vy = vy + 0.5f * (ay - x * vy) * dt;
		vz = vz + 0.5f * (az - x * vz) * dt;
		vel[idx * 4 + 0] = vx;
		vel[idx * 4 + 1] = vy;
		vel[idx * 4 + 2] = vz;
	}
}

// RawNVT fusion SAFE fused kernel: BIT-IDENTICAL to (nvt_finalize_xi ; nvt_second_step).
// stage 1 reduces the KE partials with the SAME 256-wide tree as nvt_finalize_xi
// (only threads 0..255 participate, identical strided sum + reduction order), so
// the reduced temperature -> xi update is bit-for-bit identical to the 4-kernel
// path.  stage 3 then applies the SAME per-member secondstep math using ALL 1024
// threads.  Launched as ONE block of _FUSE_THREADS so __syncthreads orders the xi
// update before every velocity update.  Verified bit-exact on GPU (xi and vel).
extern "C" __global__ void nvt_finalize_xi_secondstep_safe(
		const float* __restrict__ partial, const int nblocks,
		float* __restrict__ xi, const float dt_over_tau2,
		const float inv_nfreedom, const float target_temp,
		const int nme, const int* __restrict__ member,
		float* __restrict__ vel, float* __restrict__ force, const float dt)
{
	__shared__ float sm[256];
	__shared__ float sh_xi;
	int tx = threadIdx.x;

	// stage 1: 256-wide reduction IDENTICAL to nvt_finalize_xi (threads 0..255).
	if (tx < 256) {
		float s = 0.0f;
		for (int o = tx; o < nblocks; o += 256) s += partial[o];
		sm[tx] = s;
	}
	__syncthreads();
	for (int offs = 128; offs > 0; offs >>= 1) {
		if (tx < offs) sm[tx] += sm[tx + offs];
		__syncthreads();
	}
	// stage 2: thread 0 applies the SAME thermostat update as nvt_finalize_xi.
	if (tx == 0) {
		float temp = sm[0] * inv_nfreedom;
		float xnew = xi[0] + dt_over_tau2 * (temp / target_temp - 1.0f);
		xi[0] = xnew;
		xi[1] = temp;       // debug cache: last device-side temperature
		sh_xi = xnew;
	}
	__syncthreads();        // new xi visible to all threads before use
	float x = sh_xi;

	// stage 3: SAME per-member secondstep math as nvt_second_step, all threads.
	for (int i = tx; i < nme; i += blockDim.x) {
		int idx = member[i];
		float m  = vel[idx * 4 + 3];
		float ax = force[idx * 3 + 0] / m;
		float ay = force[idx * 3 + 1] / m;
		float az = force[idx * 3 + 2] / m;
		force[idx * 3 + 0] = ax;
		force[idx * 3 + 1] = ay;
		force[idx * 3 + 2] = az;
		float vx = vel[idx * 4 + 0];
		float vy = vel[idx * 4 + 1];
		float vz = vel[idx * 4 + 2];
		vx = vx + 0.5f * (ax - x * vx) * dt;
		vy = vy + 0.5f * (ay - x * vy) * dt;
		vz = vz + 0.5f * (az - x * vz) * dt;
		vel[idx * 4 + 0] = vx;
		vel[idx * 4 + 1] = vy;
		vel[idx * 4 + 2] = vz;
	}
}
'''


def _import_cupy():
	'''Lazy import of cupy via the central CuPy bridge (clear error message).'''
	return import_cupy("RawKernel NVT")


def _profile_enabled():
	v = os.environ.get("PYGAMD_RAWNVT_PROFILE", "")
	return v not in ("", "0", "false", "False", "no", "off")


def _numba_compat_enabled():
	v = os.environ.get("PYGAMD_RAWNVT_NUMBA_COMPAT", "")
	return v not in ("", "0", "false", "False", "no", "off", "NO", "OFF")


# --- RawNVT fusion RawNVT fusion: env + sentinel force-enable --------------------
# Like the engine timing / launch-trace sentinels, PYGAMD_RAWNVT_FUSION may not
# propagate to the final CuPy subprocess (external profiler env passing is
# unreliable).  The RawNVT fusion runner therefore also drops "force RawNVT fusion"
# SENTINEL files on disk; if any is present the SAFE fused path is used.  Resolved
# at first real use (in _ensure_init), so the default path pays nothing per step.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))   # .../pygamd/integrations
_PKG_DIR = os.path.dirname(_THIS_DIR)                      # .../pygamd
_ROOT = os.path.dirname(_PKG_DIR)                          # project root

_FUSION_SAFE_VALUES = ("1", "true", "True", "TRUE", "on", "On", "ON",
                       "yes", "Yes", "YES", "safe", "Safe", "SAFE")
_FUSION_AGGR_VALUES = ("2", "aggressive", "Aggressive", "AGGRESSIVE")
_FUSION_OFF_VALUES = ("", "0", "false", "False", "no", "off", "No", "Off", "NO", "OFF")


def _fusion_flag_paths():
	paths = [
		os.path.join(_PKG_DIR, "_pygamd_force_rawnvt_fusion.flag"),
		os.path.join(_ROOT, ".pygamd_force_rawnvt_fusion.flag"),
		"/tmp/pygamd_force_rawnvt_fusion.flag",
	]
	try:
		paths.append("/tmp/pygamd_force_rawnvt_fusion_%d.flag" % os.getuid())
	except Exception:
		pass
	return paths


def _read_fusion_sentinel():
	'''True if any RawNVT-fusion sentinel file exists (cwd-independent).'''
	for p in _fusion_flag_paths():
		try:
			if os.path.exists(p):
				return True
		except Exception:
			pass
	return False


def _fusion_force_flag_active():
	return _read_fusion_sentinel()


def _resolve_fusion_mode():
	'''Return (mode, note).  mode in {"off","safe","aggressive"}.

	Priority: an explicit valid env value wins; otherwise an OFF/unset env with a
	sentinel present -> "safe"; an ILLEGAL non-empty env value -> "off" (fusion NOT
	enabled) with a diagnostic note, unless a sentinel forces "safe".  Default
	(unset, no sentinel) -> "off".'''
	raw = os.environ.get("PYGAMD_RAWNVT_FUSION", "")
	v = raw.strip() if isinstance(raw, str) else raw
	if v in _FUSION_AGGR_VALUES:
		return "aggressive", None
	if v in _FUSION_SAFE_VALUES:
		return "safe", None
	if v in _FUSION_OFF_VALUES:
		if _read_fusion_sentinel():
			return "safe", None
		return "off", None
	# illegal non-empty value: do NOT silently enable fusion
	if _read_fusion_sentinel():
		return ("safe",
		        "PYGAMD_RAWNVT_FUSION=%r is not a valid value; sentinel present -> safe" % raw)
	return ("off",
	        "PYGAMD_RAWNVT_FUSION=%r is not a valid value; falling back to off "
	        "(RawNVT fusion NOT enabled). Valid: 1/true/on/yes/safe, 2/aggressive, 0/off." % raw)


def _fusion_mode():
	return _resolve_fusion_mode()[0]


class nose_hoover:
	'''RawKernel-native Nose-Hoover NVT (same constructor signature as the
	Numba nose_hoover, so integration.nvt can swap it in transparently).'''

	_PROFILE_PERIOD = 1000

	def __init__(self, info, ci, tau, temp):
		self.info = info
		# CuPy-native device-data layer (MVP-1): core particle CuPy arrays come
		# from CuInfo (info.cu) when present, instead of scattered Numba wraps.
		self.cu = getattr(info, "cu", None)
		self.ci = ci                 # kept for dump.data; member arr is transitional
		self.tau = tau
		self.temp = temp
		self.block_size = _NVT_BLOCK
		self.xi = 0.0                # host debug cache only; device d_xi is main

		self.cp = None
		self.module = None
		self.k_first = None
		self.k_ke = None
		self.k_finalize = None
		self.k_finalize_tptr = None
		self.k_second = None
		self.k_fused = None              # aggressive fused finalize_xi + secondstep
		self.k_fused_safe = None         # RawNVT fusion bit-identical fused finalize_xi + secondstep
		self.fusion_mode, self._fusion_note = _resolve_fusion_mode()  # re-resolved in _ensure_init
		self.d_xi = None
		self.d_partial = None
		self.d_graph_target_temp = None
		self.nme = int(ci.ps.nme)
		self.nfreedom = int(ci.ps.nfreedom)
		self.nblocks_me = math.ceil(self.nme / self.block_size)
		self._wrap_cache = {}

		# profiler
		self._profile = _profile_enabled()
		self._calls = 0
		self._temp_reduce = 0
		self._xi_updates = 0
		self._last_target = 0.0
		self.numba_compat = _numba_compat_enabled()

	# --- user API -------------------------------------------------------------
	def setT(self, temp):
		self.temp = temp

	def _get_target_temp(self, timestep):
		if hasattr(self.temp, "getValue"):
			return float(self.temp.getValue(timestep))
		return float(self.temp)

	def register(self, timestep):
		# Default RawKernel NVT keeps temperature reduction on device.  The opt-in
		# numba compatibility path asks chare.comp_info for the same host-side
		# temperature/cache behavior as integrations/nvt.py in pygamd-numba.
		if self.numba_compat:
			self.info.compute_properties['temperature'] = True

	# --- lazy device setup ----------------------------------------------------
	def _ensure_init(self):
		if self.cp is not None:
			return
		# RawNVT fusion: re-resolve the fusion mode from env + sentinel at first real
		# use.  The runner's sentinel may be created after this object was built but
		# before the timed run; this runs ONCE (cp is set below), so the default
		# path pays no per-step cost.
		self.fusion_mode, self._fusion_note = _resolve_fusion_mode()
		cp = _import_cupy()
		self.cp = cp
		try:
			mod = cp.RawModule(code=_KERNEL_SRC)
			self.k_first = mod.get_function("nvt_first_step")
			self.k_ke = mod.get_function("nvt_ke_partial")
			self.k_finalize = mod.get_function("nvt_finalize_xi")
			self.k_finalize_tptr = mod.get_function("nvt_finalize_xi_tptr")
			self.k_second = mod.get_function("nvt_second_step")
			self.k_fused = mod.get_function("nvt_finalize_xi_secondstep")
			self.k_fused_safe = mod.get_function("nvt_finalize_xi_secondstep_safe")
			self.module = mod
		except Exception as e:
			raise RuntimeError(
				"Failed to compile RawKernel NVT kernels. Original error: %s" % repr(e)
			) from e
		# device thermostat scalar: [0]=xi, [1]=last temperature (debug)
		self.d_xi = cp.zeros(2, dtype=cp.float32)
		self.d_partial = cp.zeros(max(self.nblocks_me, 1), dtype=cp.float32)
		self.d_graph_target_temp = cp.zeros((1,), dtype=cp.float32)

	def _wrap(self, key, arr):
		cached = self._wrap_cache.get(key)
		if cached is not None and cached[0] is arr:
			return cached[1]
		w = as_cupy_array(arr, ravel=True)
		self._wrap_cache[key] = (arr, w)
		return w

	# --- CuInfo-preferred core-array accessor ---------------------------------
	def _dev(self, name, fallback_arr):
		'''1-D CuPy view of a core particle array, preferring the CuPy-native
		CuInfo device-data layer (centralized); transitional path wraps the Numba
		array zero-copy.  Both yield a CuPy array for the same RawKernel.'''
		if self.cu is not None:
			return self.cu.dev(name)
		return self._wrap(name, fallback_arr)

	# --- integration steps ----------------------------------------------------
	def firststep(self, timestep):
		self._ensure_init()
		if _sp.ENABLED: _sp.inc("integration_firststep_count")
		member = self._wrap("member", self.ci.ps.d_member)
		pos    = self._dev("pos", self.info.d_pos)
		vel    = self._dev("vel", self.info.d_vel)
		force  = self._dev("force", self.info.d_force)
		image  = self._dev("image", self.info.d_image)
		box    = self._dev("box", self.info.d_box)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks_me),)
		block = (int(self.block_size),)
		prof.begin("rawnvt_firststep")
		prof.begin("rawnvt_total")
		if _lt.LAUNCH: _lt.launch("integrator_nvt_firststep")
		if _t.ENABLED: _t.gpu_begin("integrator_nvt_firststep")
		self.k_first(grid, block, (np.int32(self.info.npa), np.int32(self.nme), member,
		                           pos, vel, force, image, box, self.d_xi, dt))
		if _t.ENABLED: _t.gpu_end("integrator_nvt_firststep")
		prof.end("rawnvt_total")
		prof.end("rawnvt_firststep")

	def secondstep(self, timestep):
		self._ensure_init()
		if _sp.ENABLED:
			_sp.inc("integration_secondstep_count")
			# RawNVT fusion: surface the (opt-in) fusion mode for the launch-attribution
			# report.  Idempotent; reflects the env+sentinel resolution done in
			# _ensure_init.
			_sp.set_meta("rawnvt_fusion_mode", self.fusion_mode)
			_sp.set_meta("rawnvt_fusion_enabled", self.fusion_mode != "off")
			_sp.set_meta("rawnvt_fusion_force_flag_active", _fusion_force_flag_active())
			if self._fusion_note:
				_sp.note(self._fusion_note)
		target_temp = self._get_target_temp(timestep)
		if target_temp <= 0.0:
			raise RuntimeError("NVT target temperature must be positive")
		self._last_target = target_temp
		if self.numba_compat:
			return self._secondstep_numba_compat(timestep, target_temp)

		member = self._wrap("member", self.ci.ps.d_member)
		vel    = self._dev("vel", self.info.d_vel)
		force  = self._dev("force", self.info.d_force)
		dt = np.float32(self.info.dt)
		dt_over_tau2 = np.float32(self.info.dt / (self.tau * self.tau))
		inv_nfreedom = np.float32(1.0 / float(self.nfreedom))
		grid = (int(self.nblocks_me),)
		block = (int(self.block_size),)

		prof.begin("rawnvt_total")
		# 1) device-side kinetic partial sums (multi-block)
		prof.begin("rawnvt_temp_partial")
		prof.count("rawnvt_temp_reduce_calls")
		if _lt.LAUNCH: _lt.launch("integrator_nvt_ke_reduce")
		if _t.ENABLED: _t.gpu_begin("integrator_nvt_ke_reduce")
		self.k_ke(grid, block, (np.int32(self.nme), member, vel, self.d_partial))
		if _t.ENABLED: _t.gpu_end("integrator_nvt_ke_reduce")
		prof.end("rawnvt_temp_partial")

		if self.fusion_mode != "off":
			# RawNVT fusion 2+3 fused: reduce -> xi update -> secondstep in ONE single
			# block (4 -> 3 launches/step).  "safe" = bit-identical kernel (default
			# when enabled); "aggressive" = legacy aggressive fused kernel (regression only).
			fused_k = self.k_fused_safe if self.fusion_mode == "safe" else self.k_fused
			prof.begin("rawnvt_xi_secondstep")
			prof.count("rawnvt_xi_update_calls")
			if _lt.LAUNCH: _lt.launch("integrator_nvt_xi_secondstep")
			if _t.ENABLED: _t.gpu_begin("integrator_nvt_xi_secondstep")
			fused_k((1,), (_FUSE_THREADS,),
			        (self.d_partial, np.int32(self.nblocks_me), self.d_xi,
			         dt_over_tau2, inv_nfreedom, np.float32(target_temp),
			         np.int32(self.nme), member, vel, force, dt))
			if _t.ENABLED: _t.gpu_end("integrator_nvt_xi_secondstep")
			prof.end("rawnvt_xi_secondstep")
		else:
			# legacy 4-kernel path (PYGAMD_RAWNVT_FUSION=0)
			prof.begin("rawnvt_xi_update")
			prof.count("rawnvt_xi_update_calls")
			if _lt.LAUNCH: _lt.launch("integrator_nvt_xi_update")
			if _t.ENABLED: _t.gpu_begin("integrator_nvt_xi_update")
			self.k_finalize((1,), (_FINALIZE_THREADS,),
			                (self.d_partial, np.int32(self.nblocks_me), self.d_xi,
			                 dt_over_tau2, inv_nfreedom, np.float32(target_temp)))
			if _t.ENABLED: _t.gpu_end("integrator_nvt_xi_update")
			prof.end("rawnvt_xi_update")
			prof.begin("rawnvt_secondstep")
			if _lt.LAUNCH: _lt.launch("integrator_nvt_secondstep")
			if _t.ENABLED: _t.gpu_begin("integrator_nvt_secondstep")
			self.k_second(grid, block, (np.int32(self.nme), member, vel, force, self.d_xi, dt))
			if _t.ENABLED: _t.gpu_end("integrator_nvt_secondstep")
			prof.end("rawnvt_secondstep")
		prof.end("rawnvt_total")

		if self._profile:
			self._calls += 1
			self._temp_reduce += 1
			self._xi_updates += 1
			if self._calls % self._PROFILE_PERIOD == 0:
				print("info : [rawnvt profile] calls=%d temp_reduce=%d xi_updates=%d target_temp=%.6g"
				      % (self._calls, self._temp_reduce, self._xi_updates, self._last_target))

	def _secondstep_numba_compat(self, timestep, target_temp):
		self.ci.calculate(timestep)
		self.xi += self.info.dt / (self.tau * self.tau) * (self.ci.temp / target_temp - 1.0)
		self.d_xi[0] = np.float32(self.xi)
		self.d_xi[1] = np.float32(self.ci.temp)
		member = self._wrap("member", self.ci.ps.d_member)
		vel    = self._dev("vel", self.info.d_vel)
		force  = self._dev("force", self.info.d_force)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks_me),)
		block = (int(self.block_size),)
		prof.begin("rawnvt_total")
		prof.begin("rawnvt_secondstep")
		if _lt.LAUNCH: _lt.launch("integrator_nvt_secondstep")
		if _t.ENABLED: _t.gpu_begin("integrator_nvt_secondstep")
		self.k_second(grid, block, (np.int32(self.nme), member, vel, force, self.d_xi, dt))
		if _t.ENABLED: _t.gpu_end("integrator_nvt_secondstep")
		prof.end("rawnvt_secondstep")
		prof.end("rawnvt_total")

	# --- debug helper (optional; performs a host sync, not used per-step) ------
	def get_temperature(self):
		self._ensure_init()
		return float(asnumpy(self.d_xi)[1])

	# --- CUDA-graph helpers (launch-only; no prof sync) ---------------
	# Graph replay captures the pointer to d_graph_target_temp.  graph_prepare_replay()
	# updates that one-element device scalar before each launch, so VariantLinear and
	# other variant targets are graph-safe instead of forcing variant_temp fallback.
	def _graph_temp_is_variant(self):
		return hasattr(self.temp, "getValue")

	def _target_for_step(self, timestep):
		target = float(self._get_target_temp(timestep))
		if target <= 0.0:
			raise RuntimeError("NVT target temperature must be positive")
		return target

	def _ensure_graph_target_temp(self):
		self._ensure_init()
		if self.d_graph_target_temp is None:
			self.d_graph_target_temp = self.cp.zeros((1,), dtype=self.cp.float32)
		return self.d_graph_target_temp

	def graph_ready(self):
		if self.numba_compat:
			return False
		try:
			self._ensure_init()
			self._ensure_graph_target_temp()
			self._target_for_step(0)
		except Exception:
			return False
		return (self.k_second is not None and self.k_first is not None
		        and self.k_finalize_tptr is not None and self.d_graph_target_temp is not None)

	def graph_unready_reason(self):
		if self.numba_compat:
			return "nvt_numba_compat_host_temperature"
		try:
			self._target_for_step(0)
		except Exception:
			return "nvt_invalid_target_temperature"
		if self.cp is None or self.k_first is None or self.k_second is None or self.k_finalize_tptr is None:
			return "nvt_not_initialized"
		return "unsupported_integrator"

	def graph_target_value(self):
		# A stable signature target only.  The captured graph reads the live value from
		# d_graph_target_temp, updated in graph_prepare_replay(timestep).
		return float(self._target_for_step(0))

	def graph_prepare_replay(self, timestep):
		d_temp = self._ensure_graph_target_temp()
		target = np.float32(self._target_for_step(timestep))
		d_temp.fill(target)
		if _sp.ENABLED:
			_sp.inc("cuda_graph_temperature_scalar_update_count")
			_sp.inc("cuda_graph_nvt_temperature_scalar_update_count")
			if self._graph_temp_is_variant():
				_sp.inc("cuda_graph_variant_temperature_update_count")
				_sp.inc("cuda_graph_nvt_variant_temperature_update_count")

	def graph_signature(self):
		mode = "variant_temperature" if self._graph_temp_is_variant() else "static_temperature"
		try:
			target0 = float(self._target_for_step(0))
		except Exception:
			target0 = "invalid"
		return ("nvt", int(self.nme), int(self.nfreedom), int(self.block_size),
		        float(self.info.dt), float(self.tau), mode, id(self.temp) if self._graph_temp_is_variant() else target0,
		        target0, id(self.d_xi), id(self.d_partial), id(self._ensure_graph_target_temp()))

	def graph_firststep(self):
		'''Launch ONLY the NVT firststep kernel (for full-graph capture/replay).
		Same kernel and semantics as firststep(); no prof, no host sync.'''
		member = self._wrap("member", self.ci.ps.d_member)
		pos    = self._dev("pos", self.info.d_pos)
		vel    = self._dev("vel", self.info.d_vel)
		force  = self._dev("force", self.info.d_force)
		image  = self._dev("image", self.info.d_image)
		box    = self._dev("box", self.info.d_box)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks_me),); block = (int(self.block_size),)
		self.k_first(grid, block, (np.int32(self.info.npa), np.int32(self.nme), member,
		                           pos, vel, force, image, box, self.d_xi, dt))

	def graph_temp_partial(self):
		member = self._wrap("member", self.ci.ps.d_member)
		vel    = self._dev("vel", self.info.d_vel)
		grid = (int(self.nblocks_me),); block = (int(self.block_size),)
		self.k_ke(grid, block, (np.int32(self.nme), member, vel, self.d_partial))

	def graph_xi_update(self, target_temp=None):
		dt_over_tau2 = np.float32(self.info.dt / (self.tau * self.tau))
		inv_nfreedom = np.float32(1.0 / float(self.nfreedom))
		self.k_finalize_tptr((1,), (_FINALIZE_THREADS,),
		                    (self.d_partial, np.int32(self.nblocks_me), self.d_xi,
		                     dt_over_tau2, inv_nfreedom, self._ensure_graph_target_temp()))

	def graph_secondstep(self):
		member = self._wrap("member", self.ci.ps.d_member)
		vel    = self._dev("vel", self.info.d_vel)
		force  = self._dev("force", self.info.d_force)
		dt = np.float32(self.info.dt)
		grid = (int(self.nblocks_me),); block = (int(self.block_size),)
		self.k_second(grid, block, (np.int32(self.nme), member, vel, force, self.d_xi, dt))
