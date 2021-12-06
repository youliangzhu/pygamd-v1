import pygamd
from numba import cuda
import math
	
mst = pygamd.snapshot.read("lj.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

fn = pygamd.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="a", type_j="a", alpha=30.0, sigma=3.0)
app.add(fn)

# inn = pygamd.integration.gwvv(info=mst, group=['a'])
# app.add(inn)

di = pygamd.dump.data(info=mst, group=['a'], file='data.log', period=100)
app.add(di)

app.run(1000)