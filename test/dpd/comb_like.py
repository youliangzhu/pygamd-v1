import pygamd
	
mst = pygamd.snapshot.read("comblike.mst")
app = pygamd.application.dynamics(info=mst, dt=0.04)

fn = pygamd.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
fn.setParams(type_i="S", type_j="S", alpha=25.0, sigma=3.0)
fn.setParams(type_i="A", type_j="S", alpha=150.0, sigma=3.0)
fn.setParams(type_i="B", type_j="S", alpha=27.0, sigma=3.0)
fn.setParams(type_i="A", type_j="B", alpha=70.0, sigma=3.0)
app.add(fn)

fb = pygamd.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type = 'A-A', param=[4.0, 0.0])#(,K, R0)
fb.setParams(bond_type = 'A-B', param=[4.0, 0.0])#(,K, R0)
fb.setParams(bond_type = 'B-B', param=[4.0, 0.0])#(,K, R0)
app.add(fb)

inn = pygamd.integration.gwvv(info=mst, group=['A', 'B', 'S'])
app.add(inn)

di = pygamd.dump.data(info=mst, group=['A', 'B', 'S'], file='data.log', period=500)
app.add(di)

dm = pygamd.dump.mst(info=mst, group=['A', 'B', 'S'], file='p.mst', period=100000)
app.add(dm)

app.run(4000000)
