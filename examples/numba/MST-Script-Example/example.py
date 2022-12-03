pygamd pygamd

mst = pygamd.snapshot.read("example.mst")
app = pygamd.application.dynamics(info=mst, dt=0.001)

fn = pygamd.force.nonbonded(info=mst, rcut=3.0, func='lj')
fn.setParams(type_i="A", type_j="A", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

fb = pygamd.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type = 'A-A', param=[1000.0, 0.5])#(,k, r0)
app.add(fb)

fa = pygamd.force.angle(info=mst, func='harmonic')
fa.setParams(angle_type='A-A-A', param=[100.0, 120.0])#(,k, t0)
app.add(fa)

 
fd = pygamd.force.dihedral(info=mst, func='harmonic')
fd.setParams(dihedral_type = 'A-A-A-A', param=[25.0, 180.0])#(,k, t0)
app.add(fd)

nvt = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=1.0)
app.add(nvt)

dd = pygamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

#run 
app.run(20000)

