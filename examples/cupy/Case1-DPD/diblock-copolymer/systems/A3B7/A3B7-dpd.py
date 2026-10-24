# CuPy-native version of the A5B5 DPD simulation.
import pygamd


mst = pygamd.snapshot.read("A3B7.mst")
app = pygamd.application.dynamics(info=mst, dt=0.04, sort=True)

fn = pygamd.force.dpd(info=mst, rcut=1.0)
fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
app.add(fn)

fb = pygamd.force.bond(info=mst, func="harmonic")
fb.setParams(bond_type="A-A", param=[4.0, 0.0])
fb.setParams(bond_type="A-B", param=[4.0, 0.0])
fb.setParams(bond_type="B-B", param=[4.0, 0.0])
app.add(fb)

gw = pygamd.integration.gwvv(info=mst, group="all", dpd_force=fn)
app.add(gw)


di = pygamd.dump.data(info=mst, group="all", file="data.log", period=500)
app.add(di)

dm = pygamd.dump.mst(info=mst, group="all", file="p.mst", period=100000)
app.add(dm)


xml = pygamd.dump.xml(info=mst, group="all", file="p", period=10000)
app.add(xml)

app.run(500000)
