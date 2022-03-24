# the script for the simulation of Lennard-Jones liquids
import pygamd
from numba import cuda
import numba as nb
	
mst = pygamd.snapshot.read("lj.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

# user-defined potential form
@cuda.jit(device=True)
def lj(rsq, param, fp):
	epsilon = param[0]
	sigma = param[1]
	alpha = param[2]
	rcut = param[3]
	if rsq<rcut*rcut:
		sigma2 = sigma*sigma
		r2inv = sigma2/rsq;
		r6inv = r2inv * r2inv * r2inv;
		f = nb.float32(4.0) * r2inv * r6inv * (nb.float32(12.0) * r6inv - nb.float32(6.0) * alpha)/sigma2	
		p = nb.float32(4.0) * r6inv * ( r6inv - nb.float32(1.0))
		fp[0]=f
		fp[1]=p	

# non-bonded interactions
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

