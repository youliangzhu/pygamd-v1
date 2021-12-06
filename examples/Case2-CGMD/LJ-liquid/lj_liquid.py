import pygamd
	
mst = pygamd.snapshot.read("lj.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

fn = pygamd.force.nonbonded(info=mst, rcut=3.0, func='lj')
fn.setParams(type_i="a", type_j="a", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

inn = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=1.0, temperature=1.0)
app.add(inn)

dd = pygamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

app.run(20000)
