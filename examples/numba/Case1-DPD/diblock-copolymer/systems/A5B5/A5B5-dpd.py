import sys
sys.path.insert(0, "/home/aswinkumar99/pygamd-v1")

import pygamd

print("Using pygamd from:", pygamd.__path__)

mst = pygamd.snapshot.read("A5B5.mst")
app = pygamd.application.dynamics(info=mst, dt=0.04)

print(type(mst))
fn = pygamd.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
app.add(fn)

fb = pygamd.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type = 'A-A', param=[4.0, 0.0])#param=[k, r0]
fb.setParams(bond_type = 'A-B', param=[4.0, 0.0])#param=[k, r0]
fb.setParams(bond_type = 'B-B', param=[4.0, 0.0])#param=[k, r0]
app.add(fb)

gw = pygamd.integration.gwvv(info=mst, group='all')
app.add(gw)


di = pygamd.dump.data(info=mst, group='all', file='data.log', period=500)
app.add(di)

dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=100000)
app.add(dm)

#ready ro run
app.run(500000)#(the number of steps to run)
