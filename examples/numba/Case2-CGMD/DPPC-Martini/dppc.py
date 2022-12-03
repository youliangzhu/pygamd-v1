import pygamd
from numba import cuda
import numba as nb
import math

mst = pygamd.snapshot.read("dppc.mst")
app = pygamd.application.dynamics(info=mst, dt=0.02)

fn = pygamd.force.nonbonded_c(info=mst, rcut=1.2, func='lj_coulomb')
fn.setParams(type_i='P4', type_j='P4', param=[5.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Na', type_j='Na', param=[4.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='C1', type_j='C1', param=[3.5, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Qa', type_j='Qa', param=[5.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Q0', type_j='Q0', param=[3.5, 0.47, 1.0, 15.0, 1.2]) 
fn.setParams(type_i='P4', type_j='Na', param=[4.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='P4', type_j='C1', param=[2.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='P4', type_j='Qa', param=[5.6, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='P4', type_j='Q0', param=[5.6, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Na', type_j='C1', param=[2.7, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Na', type_j='Qa', param=[4.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='Na', type_j='Q0', param=[4.0, 0.47, 1.0, 15.0, 1.2])
fn.setParams(type_i='C1', type_j='Qa', param=[2.0, 0.62, 1.0, 15.0, 1.2]) 
fn.setParams(type_i='C1', type_j='Q0', param=[2.0, 0.62, 1.0, 15.0, 1.2])
fn.setParams(type_i='Qa', type_j='Q0', param=[4.5, 0.47, 1.0, 15.0, 1.2])
app.add(fn)

fb = pygamd.force.bond(info=mst, func="harmonic")
fb.setParams(bond_type = 'Q0-Qa', param=[1250.0, 0.47])#param=[k, r0]
fb.setParams(bond_type = 'Qa-Na', param=[1250.0, 0.47])#param=[k, r0]
fb.setParams(bond_type = 'Na-Na', param=[1250.0, 0.37])#param=[k, r0]
fb.setParams(bond_type = 'Na-C1', param=[1250.0, 0.47])#param=[k, r0]
fb.setParams(bond_type = 'C1-C1', param=[1250.0, 0.47])#param=[k, r0]
app.add(fb)

fa = pygamd.force.angle(info=mst, func="harmonic_cos")
fa.setParams(angle_type='Qa-Na-Na', param=[25.0, 120.0])#param=[k, t0]
fa.setParams(angle_type='Qa-Na-C1', param=[25.0, 180.0])#param=[k, t0]
fa.setParams(angle_type='Na-C1-C1', param=[25.0, 180.0])#param=[k, t0]
fa.setParams(angle_type='C1-C1-C1', param=[25.0, 180.0])#param=[k, t0]
app.add(fa)

Temp = 318.000  #k
Temp*=8.3143/1000.0#reduced unit
nvt = pygamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=Temp)
app.add(nvt)

dd = pygamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

xml = pygamd.dump.xml(info=mst, group='all', file='p', period=10000)
app.add(xml)

#run 
app.run(2000000)

