import pygamd

mst = pygamd.snapshot.read("A1B23.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

slj = pygamd.force.slj(info=mst, rcut=3.0)
slj.setParams(type_i="A", type_j="A", params=[1.0,1.0,1.0]) # [epsilon, alpha, sigma, rcut]
slj.setParams(type_i="A", type_j="B", params=[1.0,1.0,1.0])
slj.setParams(type_i="B", type_j="B", params=[1.0,1.0,1.0])
app.add(slj)

fb = pygamd.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type='A-B', param=[1000.0, 2.0])  # param=[k, r0]
fb.setParams(bond_type='B-B', param=[1000.0, 1.0])  # param=[k, r0]
app.add(fb)

inn = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=1.0, temperature=1.0)
app.add(inn)

di = pygamd.dump.data(info=mst, group='all', file='data.log', period=500)
app.add(di)

dm = pygamd.dump.xml(info=mst, group='all', file='p.xml', period=500, properties=['position', 'type', 'bond', 'diameter'])
app.add(dm)

# ready ro run
app.run(50000)  # (the number of steps to run)
