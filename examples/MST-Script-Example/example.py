import jugamd

mst = jugamd.snapshot.read("example.mst")
app = jugamd.application.dynamics(info=mst, dt=0.001)

fn = jugamd.force.nonbonded(info=mst, rcut=3.0, func='lj')
fn.setParams(type_i="A", type_j="A", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

fb = jugamd.force.bond(info=mst, func='harmonic')
fb.setParams(bond_type = 'A-A', param=[1000.0, 0.5])#(,k, r0)
app.add(fb)

fa = jugamd.force.angle(info=mst, func='harmonic')
fa.setParams(angle_type='A-A-A', param=[100.0, 120.0])#(,k, t0)
app.add(fa)

 
fd = jugamd.force.dihedral(info=mst, func='harmonic')
fd.setParams(dihedral_type = 'A-A-A-A', param=[25.0, 180.0])#(,k, t0)
app.add(fd)

nvt = jugamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=1.0)
app.add(nvt)

dd = jugamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = jugamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

#run 
app.run(20000)

