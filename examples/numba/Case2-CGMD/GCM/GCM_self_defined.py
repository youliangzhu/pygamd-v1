import pygamd
from numba import cuda
import numba as nb
import numpy as np
import math

mst = pygamd.snapshot.read("init_one.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

@cuda.jit(device=True)
def gem(rsq, param, fp):
	epsilon = param[0]
	sigma = param[1]
	n = param[2]
	rcut = param[3]
	if rsq<rcut*rcut:
		r = math.sqrt(rsq)
		val_p = math.pow(r/sigma, n)
		val_e = math.exp(-val_p)

		f = np.float32(0.0)
		if (rsq>np.float32(0.0)):
			f= epsilon*val_e*val_p*n/rsq
		p = epsilon*val_e
		fp[0]=f
		fp[1]=p	

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