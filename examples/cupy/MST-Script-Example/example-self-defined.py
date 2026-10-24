import pygamd
from numba import cuda
import numba as nb
import math


mst = pygamd.snapshot.read("example.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

LJ_CODE = r"""
__device__ __forceinline__ void lj_self(float rsq, const float* param, float* fp) {
    float epsilon = param[0];
    float sigma = param[1];
    float alpha = param[2];
    float rcut = param[3];
    if (rsq < rcut * rcut) {
        float sigma2 = sigma * sigma;
        float r2inv = sigma2 / rsq;
        float r6inv = r2inv * r2inv * r2inv;
        float f = 4.0f * epsilon * r2inv * r6inv * (12.0f * r6inv - 6.0f * alpha) / sigma2;
        float p = 4.0f * epsilon * r6inv * (r6inv - 1.0f);
        fp[0] = f;
        fp[1] = p;
    }
}
"""

lj = pygamd.potentials.pair_raw(name="lj_self", code=LJ_CODE, nparams=4)
fn = pygamd.force.nonbonded(info=mst, rcut=3.0, func=lj)
fn.setParams(type_i="A", type_j="A", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

BOND_CODE = r"""
__device__ __forceinline__ void bond_harmonic_self(float rsq, const float* param, float* fp) {
    float k = param[0];
    float r0 = param[1];
    float r = sqrtf(rsq);
    float f = k * (r0 / r - 1.0f);
    float p = 0.5f * k * (r0 - r) * (r0 - r);
    fp[0] = f;
    fp[1] = p;
}
"""

bond_harmonic = pygamd.potentials.bond_raw(name="bond_harmonic_self", code=BOND_CODE, nparams=2)
fb = pygamd.force.bond(info=mst, func=bond_harmonic)
fb.setParams(bond_type = 'A-A', param=[1000.0, 0.5])#param=[k, r0]
app.add(fb)


ANGLE_CODE = r"""
__device__ __forceinline__ void angle_harmonic_self(float cos_abc, float sin_abc,
                                                    const float* param, float* fp) {
    float k = param[0];
    float t0 = param[1];
    float dth = acosf(cos_abc) - 3.14159265358979323846f * t0 / 180.0f;
    float f = k * dth;
    float p = 0.5f * f * dth;
    fp[0] = f;
    fp[1] = p;
}
"""

angle_harmonic = pygamd.potentials.angle_raw(name="angle_harmonic_self", code=ANGLE_CODE, nparams=2)
fa = pygamd.force.angle(info=mst, func=angle_harmonic)
fa.setParams(angle_type='A-A-A', param=[100.0, 120.0])#param=[k, t0]
app.add(fa)

DIHEDRAL_CODE = r"""
__device__ __forceinline__ void dihedral_harmonic_self(float cos_abcd, float sin_abcd,
                                                       const float* param, float* fp) {
    float k = param[0];
    float cos_phi0 = param[1];
    float sin_phi0 = param[2];
    float cos_factor = param[3];
    float f = cos_factor * (-sin_abcd * cos_phi0 + cos_abcd * sin_phi0);
    float p = 1.0f + cos_factor * (cos_abcd * cos_phi0 + sin_abcd * sin_phi0);
    fp[0] = -k * f;
    fp[1] = k * p;
}
"""
	
dihedral_harmonic = pygamd.potentials.dihedral_raw(name="dihedral_harmonic_self", code=DIHEDRAL_CODE, nparams=4)
fd = pygamd.force.dihedral(info=mst, func=dihedral_harmonic)
fd.setParams(dihedral_type = 'A-A-A-A', param=[25.0, math.cos(math.pi), math.sin(math.pi), -1.0])#param=[k, cos_phi0, sin_phi0, cos_factor]
app.add(fd)

nvt = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=1.0)
app.add(nvt)

dd = pygamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

#run 
app.run(20000)

