import pygamd
from numba import cuda
import numba as nb

	
mst = pygamd.snapshot.read("lj.mst")



@cuda.jit(device=True)
def harmonic(rsq, param, fp):
	k = param[0]
	r0 = param[1]
	r = math.sqrt(rsq)
	f = k * (r0/r - nb.float32(1.0))
	p = nb.float32(0.5) * k * (r0 - r) * (r0 - r)
	fp[0]=f
	fp[1]=p


fn = pygamd.force.bond(info=mst, func=harmonic)
fn.setParams('a-a', [100.0, 1.0])
fn.compute(0)


