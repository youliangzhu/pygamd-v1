"""Built-in CuPy RawKernel non-bonded potential descriptors."""

from pygamd.forces.potentials.base import pair_raw, _device_potential


_LJ_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, const float* param, float* fp) {
    float epsilon = param[0];
    float sigma   = param[1];
    float alpha   = param[2];
    float rcut    = param[3];
    if (rsq < rcut * rcut) {
        float sigma2 = sigma * sigma;
        float r2inv  = sigma2 / rsq;
        float r6inv  = r2inv * r2inv * r2inv;
        float f = 4.0f * epsilon * r2inv * r6inv * (12.0f * r6inv - 6.0f * alpha) / sigma2;
        float p = 4.0f * epsilon * r6inv * (r6inv - 1.0f);
        fp[0] = f;
        fp[1] = p;
    }
}
'''


_GEM4_RAW_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, const float* param, float* fp) {
    float epsilon = param[0];
    float sigma   = param[1];
    float rcut    = param[3];
    if (rsq < rcut * rcut) {
        float inv_sigma2 = 1.0f / (sigma * sigma);
        float q     = rsq * inv_sigma2;
        float val_p = q * q;
        float val_e = expf(-val_p);
        float f = 0.0f;
        if (rsq > 0.0f) {
            f = epsilon * val_e * val_p * 4.0f / rsq;
        }
        fp[0] = f;
        fp[1] = epsilon * val_e;
    }
}
'''


_LJ_COULOMB_CHARGE_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, float qi, float qj,
        const float* param, float* fp)
{
    float lj1 = param[0];
    float lj2 = param[1];
    float coulomb_eff = param[2];
    float rcutsq = param[3];
    if (rsq < rcutsq)
    {
        float r2inv = 1.0f / rsq;
        float rinv = sqrtf(r2inv);
        float r6inv = r2inv * r2inv * r2inv;
        float qiqj = qi * qj;
        fp[0] = r2inv * r6inv * (12.0f * lj1 * r6inv - 6.0f * lj2)
              + coulomb_eff * qiqj * r2inv * rinv;
        fp[1] = r6inv * (lj1 * r6inv - lj2) + coulomb_eff * qiqj * rinv;
    }
}
'''


_LJ_COULOMB_ERFC_CHARGE_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, float qi, float qj,
        const float* param, float* fp)
{
    float lj1 = param[0];
    float lj2 = param[1];
    float coulomb_eff = param[2];
    float kappa = param[3];
    float rcutsq = param[4];
    if (rsq < rcutsq && rsq > 0.0f)
    {
        const float inv_sqrt_pi = 0.5641895835477563f;
        float r2inv = 1.0f / rsq;
        float rinv = sqrtf(r2inv);
        float r = rsq * rinv;
        float r6inv = r2inv * r2inv * r2inv;
        float qiqj = qi * qj;
        float kr = kappa * r;
        float erfc_kr = erfcf(kr);
        float exp_kr2 = expf(-kr * kr);
        float coul_f = coulomb_eff * qiqj *
            (erfc_kr * r2inv * rinv + 2.0f * kappa * inv_sqrt_pi * exp_kr2 * r2inv);
        float coul_u = coulomb_eff * qiqj * erfc_kr * rinv;
        fp[0] = r2inv * r6inv * (12.0f * lj1 * r6inv - 6.0f * lj2) + coul_f;
        fp[1] = r6inv * (lj1 * r6inv - lj2) + coul_u;
    }
}
'''


_COULOMB_ERFC_CHARGE_CODE = r'''
__device__ __forceinline__ void %(name)s(float rsq, float qi, float qj,
        const float* param, float* fp)
{
    float coulomb_eff = param[0];
    float kappa = param[1];
    float rcutsq = param[2];
    if (rsq < rcutsq && rsq > 0.0f)
    {
        const float inv_sqrt_pi = 0.5641895835477563f;
        float r2inv = 1.0f / rsq;
        float rinv = sqrtf(r2inv);
        float r = rsq * rinv;
        float qiqj = qi * qj;
        float kr = kappa * r;
        float erfc_kr = erfcf(kr);
        float exp_kr2 = expf(-kr * kr);
        fp[0] = coulomb_eff * qiqj *
            (erfc_kr * r2inv * rinv + 2.0f * kappa * inv_sqrt_pi * exp_kr2 * r2inv);
        fp[1] = coulomb_eff * qiqj * erfc_kr * rinv;
    }
}
'''


def lj_raw(name="lj", fast_math=False):
	"""Built-in Lennard-Jones pair potential.

	setParams order: [epsilon, sigma, alpha, rcut].
	"""
	return pair_raw(name=name, code=_LJ_RAW_CODE % {"name": name}, nparams=4,
	                fast_math=fast_math)


def gem4_raw(name="gem4", fast_math=False):
	"""Built-in GEM-4 pair potential.

	setParams order: [epsilon, sigma, n, rcut].  n is fixed to 4 internally;
	param[2] is kept for layout compatibility.
	"""
	return pair_raw(name=name, code=_GEM4_RAW_CODE % {"name": name}, nparams=4,
	                fast_math=fast_math)


def lj_coulomb(name="lj_coulomb", fast_math=False):
	"""Built-in charged LJ + Coulomb pair potential.

	setParams order follows pygamd-numba:
	    [epsilon, sigma, alpha, epsilonr, rcut]
	The RawKernel receives converted device params:
	    [lj1, lj2, 138.935/epsilonr, rcut^2]
	"""
	return _device_potential("pair_charge", name,
	                         _LJ_COULOMB_CHARGE_CODE % {"name": name},
	                         nparams=5, fast_math=fast_math,
	                         param_style="lj_coulomb", device_nparams=4)


def lj_coulomb_erfc(name="lj_coulomb_erfc", fast_math=False):
	"""Built-in charged LJ + Ewald erfc real-space pair potential.

	Use this with pygamd.force.pppm(...) for the reciprocal-space term.
	setParams order:
	    [epsilon, sigma, alpha, epsilonr, rcut, kappa]
	The RawKernel receives converted device params:
	    [lj1, lj2, 138.935/epsilonr, kappa, rcut^2]
	"""
	return _device_potential("pair_charge", name,
	                         _LJ_COULOMB_ERFC_CHARGE_CODE % {"name": name},
	                         nparams=6, fast_math=fast_math,
	                         param_style="lj_coulomb_erfc", device_nparams=5)


def coulomb_erfc(name="coulomb_erfc", fast_math=False):
	"""Built-in Ewald erfc real-space Coulomb pair potential.

	setParams order:
	    [epsilonr, rcut, kappa]
	The RawKernel receives converted device params:
	    [138.935/epsilonr, kappa, rcut^2]
	"""
	return _device_potential("pair_charge", name,
	                         _COULOMB_ERFC_CHARGE_CODE % {"name": name},
	                         nparams=3, fast_math=fast_math,
	                         param_style="coulomb_erfc", device_nparams=3)


def get_pair_charge_builtin(name):
	if name == "lj_coulomb":
		return lj_coulomb()
	if name == "lj_coulomb_erfc":
		return lj_coulomb_erfc()
	if name == "coulomb_erfc":
		return coulomb_erfc()
	raise RuntimeError("Unknown CuPy charged pair potential %r." % (name,))
