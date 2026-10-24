# the script for the simulation of Lennard-Jones liquids
import pygamd
from numba import cuda
import numba as nb
	
mst = pygamd.snapshot.read("lj.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

# user-defined potential form
LJ_CODE = r"""
__device__ __forceinline__ void lj_self(float rsq, const float* param, float* fp) {
    float sigma = param[1];
    float alpha = param[2];
    float rcut = param[3];
    if (rsq < rcut * rcut) {
        float sigma2 = sigma * sigma;
        float r2inv = sigma2 / rsq;
        float r6inv = r2inv * r2inv * r2inv;
        float f = 4.0f * r2inv * r6inv * (12.0f * r6inv - 6.0f * alpha) / sigma2;
        float p = 4.0f * r6inv * (r6inv - 1.0f);
        fp[0] = f;
        fp[1] = p;
    }
}
"""

# non-bonded interactions
lj = pygamd.potentials.pair_raw(name="lj_self", code=LJ_CODE, nparams=4)
fn = pygamd.force.nonbonded(info=mst, rcut=3.0, func=lj)
fn.setParams(type_i="a", type_j="a", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

# integration
inn = pygamd.integration.nvt(info=mst, group=['a'], method="nh", tau=1.0, temperature=1.0)
app.add(inn)

# the record of temperature, pressure and momentum
di = pygamd.dump.data(info=mst, group=['a'], file='data.log', period=100)
app.add(di)

# the number of timesteps
app.run(1000)

