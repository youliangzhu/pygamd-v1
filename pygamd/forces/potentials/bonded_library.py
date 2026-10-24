"""Built-in CuPy RawKernel bonded potential descriptors."""

from pygamd.forces.potentials.base import _device_potential


_BOND_HARMONIC_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, const float* param, float* fp)
{
    float k = param[0];
    float r0 = param[1];
    float r = sqrtf(rsq);
    if (r < 1e-6f) return;
    float dr = r0 - r;
    fp[0] = k * (r0 / r - 1.0f);
    fp[1] = 0.5f * k * dr * dr;
}
'''


_ANGLE_HARMONIC_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float cos_abc, float sin_abc,
        const float* param, float* fp)
{
    float k = param[0];
    float t0 = param[1];
    float dth = acosf(cos_abc) - t0;
    float f = k * dth;
    fp[0] = f;
    fp[1] = 0.5f * f * dth;
}
'''


_ANGLE_HARMONIC_COS_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float cos_abc, float sin_abc,
        const float* param, float* fp)
{
    float k = param[0];
    float t0 = param[1];
    float dth = -sin_abc * (cos_abc - t0);
    float f = k * dth;
    fp[0] = f;
    fp[1] = 0.5f * f * dth;
}
'''


_DIHEDRAL_PERIODIC_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float cos_abcd, float sin_abcd,
        const float* param, float* fp)
{
    float k = param[0];
    float cos_phi0 = param[1];
    float sin_phi0 = param[2];
    float cos_factor = param[3];
    float f = cos_factor * (-sin_abcd * cos_phi0 + cos_abcd * sin_phi0);
    float p = 1.0f + cos_factor * (cos_abcd * cos_phi0 + sin_abcd * sin_phi0);
    fp[0] = -k * f;
    fp[1] = k * p;
}
'''


_DIHEDRAL_HARMONIC_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float cos_abcd, float sin_abcd,
        const float* param, float* fp)
{
    float term = param[5];
    if (term < 0.5f)
    {
        float k = param[0];
        float cos_phi0 = param[2];
        float sin_phi0 = param[3];
        float cos_factor = param[4];
        float f = cos_factor * (-sin_abcd * cos_phi0 + cos_abcd * sin_phi0);
        float p = 1.0f + cos_factor * (cos_abcd * cos_phi0 + sin_abcd * sin_phi0);
        fp[0] = -k * f;
        fp[1] = k * p;
    }
    else
    {
        float k = param[0];
        float t0 = param[1];
        float dth = acosf(cos_abcd) - t0;
        float f = k * dth;
        fp[0] = f;
        fp[1] = 0.5f * f * dth;
    }
}
'''


def bond_harmonic_raw(name="bond_harmonic", fast_math=False):
	"""Harmonic bond: U = 0.5*k*(r-r0)^2, setParams [k, r0]."""
	return _device_potential("bond_raw", name,
	                         _BOND_HARMONIC_RAW_CODE % {"name": name},
	                         nparams=2, fast_math=fast_math,
	                         param_style="bond_harmonic")


def angle_harmonic_raw(name="angle_harmonic", fast_math=False):
	"""Harmonic angle: U = 0.5*k*(theta-theta0)^2.

	setParams order: [k, theta0_deg].
	"""
	return _device_potential("angle_raw", name,
	                         _ANGLE_HARMONIC_RAW_CODE % {"name": name},
	                         nparams=2, fast_math=fast_math,
	                         param_style="angle_harmonic")


def angle_harmonic_cos_raw(name="angle_harmonic_cos", fast_math=False):
	"""Cosine-space harmonic angle.

	setParams order: [k, theta0_deg]; theta0 is converted to cos(theta0).
	"""
	return _device_potential("angle_raw", name,
	                         _ANGLE_HARMONIC_COS_RAW_CODE % {"name": name},
	                         nparams=2, fast_math=fast_math,
	                         param_style="angle_harmonic_cos")


def dihedral_periodic_raw(name="dihedral_periodic", fast_math=False):
	"""Proper periodic dihedral: U = k*(1 - cos(phi - phi0)).

	setParams order: [k, phi0_deg], term="proper".
	"""
	return _device_potential("dihedral_raw", name,
	                         _DIHEDRAL_PERIODIC_RAW_CODE % {"name": name},
	                         nparams=2, fast_math=fast_math,
	                         param_style="dihedral_periodic",
	                         proper_only=True, device_nparams=5)


def dihedral_harmonic_raw(name="dihedral_harmonic", fast_math=False):
	"""Legacy pygamd-numba harmonic dihedral.

	setParams order: [k, t0_deg]; supports term="proper" and term="improper".
	"""
	return _device_potential("dihedral_raw", name,
	                         _DIHEDRAL_HARMONIC_RAW_CODE % {"name": name},
	                         nparams=2, fast_math=fast_math,
	                         param_style="dihedral_harmonic",
	                         device_nparams=6)


def get_bond_builtin(name):
	if name == "harmonic":
		return bond_harmonic_raw()
	raise RuntimeError("Unknown CuPy bond potential %r." % (name,))


def get_angle_builtin(name):
	if name == "harmonic":
		return angle_harmonic_raw()
	if name == "harmonic_cos":
		return angle_harmonic_cos_raw()
	raise RuntimeError("Unknown CuPy angle potential %r." % (name,))


def get_dihedral_builtin(name):
	if name == "periodic":
		return dihedral_periodic_raw()
	if name == "harmonic":
		return dihedral_harmonic_raw()
	raise RuntimeError("Unknown CuPy dihedral potential %r." % (name,))

