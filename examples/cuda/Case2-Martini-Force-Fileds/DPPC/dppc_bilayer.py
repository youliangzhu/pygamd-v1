#!/usr/bin/python
from poetry import force_field_gala
from poetry import cu_gala as gala 
from poetry import _options
 
filename = 'dppc_bilayer.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)
 
dt = 0.02
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 1.2, 0.12)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()
epsilon_r=15.0

lj = force_field_gala.LJCoulombShiftForce(all_info, neighbor_list, 1.2, 0.9, epsilon_r, "dppc_bilayer.force_field")
app.add(lj)

bondforce = force_field_gala.BondForceHarmonic(all_info, "dppc_bilayer.force_field")
app.add(bondforce)

angleforce = force_field_gala.AngleForceHarmonicCos(all_info, "dppc_bilayer.force_field")
app.add(angleforce)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

Temperature = 318.000  #k
T = Temperature*8.3143/1000.0#reduced unit
thermo = gala.NoseHooverNVT(all_info, group, comp_info, T, 0.5)
app.add(thermo)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(400)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(100000)# (period)
app.add(ZeroMomentum)
 
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(400)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(100000)# (period)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(100000)# (period)
xml.setOutputBond(True)
xml.setOutputAngle(True)
xml.setOutputVelocity(True)
xml.setOutputMass(True)
app.add(xml)

#ready ro run

app.run(  1000000)#(How many steps to run)

