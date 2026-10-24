import pygamd


mst = pygamd.snapshot.read("init_one.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

GEM_CODE = r"""
__device__ __forceinline__ void gem(float rsq, const float* param, float* fp) {
    float epsilon = param[0];
    float sigma = param[1];
    float n = param[2];
    float rcut = param[3];
    if (rsq < rcut * rcut) {
        float r = sqrtf(rsq);
        float val_p = powf(r / sigma, n);
        float val_e = expf(-val_p);

        float f = 0.0f;
        if (rsq > 0.0f) {
            f = epsilon * val_e * val_p * n / rsq;
        }
        float p = epsilon * val_e;
        fp[0] = f;
        fp[1] = p;
    }
}
"""

gem = pygamd.potentials.pair_raw(name="gem", code=GEM_CODE, nparams=4)
fn = pygamd.force.nonbonded(info=mst, rcut=2.5, func=gem)
fn.setParams(type_i="A", type_j="A", param=[1.0, 1.0, 4.0, 2.5])   # epsilon, sigma, n, rcut
app.add(fn)

inn = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=0.01)
app.add(inn)

di = pygamd.dump.data(info=mst, group='all', file='data.log', period=1000)
app.add(di)

xml = pygamd.dump.xml(info=mst, group='all', file='p', period=10000)
app.add(xml)

app.run(100000)