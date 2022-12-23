#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options
 
filename = 'ps0.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig( int(_options.gpu))
all_info = gala.AllInfo(build_method,perform_config)
 
dt = 0.04
app = gala.Application(all_info, dt)
randomnum = 12345

neighbor_list = gala.NeighborList(all_info, 1.0, 0.1)#(,rcut,rbuffer)
neighbor_list.setRCutPair("P", "P", 1.0)
neighbor_list.setRCutPair("W", "W", 1.0)
neighbor_list.setRCutPair("P", "W", 1.0)
neighbor_list.setRCutPair("C", "C", 3.64)
neighbor_list.setRCutPair("A", "A", 3.64)
neighbor_list.setRCutPair("A", "C", 3.64)
neighbor_list.setRCutPair("A", "P", 1.0)
neighbor_list.setRCutPair("A", "W", 1.0)
neighbor_list.setRCutPair("C", "P", 1.0)
neighbor_list.setRCutPair("C", "W", 1.0)

dpd = gala.DPDForce(all_info, neighbor_list, 1.0, randomnum)#(,,rcut,seed)
dpd.setParams("P", "P", 78.3, 3.0)
dpd.setParams("W", "W", 78.3, 3.0)
dpd.setParams("C", "C", 78.3, 3.0)
dpd.setParams("A", "A", 78.3, 3.0)
dpd.setParams("A", "C", 78.3, 3.0)
dpd.setParams("A", "P", 78.3, 3.0)
dpd.setParams("A", "W", 78.3, 3.0)
dpd.setParams("C", "P", 78.3, 3.0)
dpd.setParams("C", "W", 78.3, 3.0)
dpd.setParams("P", "W", 80.42, 3.0)
app.add(dpd)
 
bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('P-P', 4.0, 0.0)#(,K0, R0)
app.add(bondforce)

group_charge = gala.ParticleSet(all_info, "charge")
kappa=0.2
enuf = gala.ENUFForce(all_info, neighbor_list, group_charge)
enuf.setParams(kappa, 2, 2, 20, 20, 20)
app.add(enuf)

ewald = gala.DPDEwaldForce(all_info, neighbor_list, group_charge, 3.64)#(,,rcut)
ewald.setParams(kappa)
app.add(ewald)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

Gwvv = gala.DPDGWVV(all_info, group)
app.add(Gwvv)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(80)
app.add(sort_method) 

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(0)# (period)
app.add(ZeroMomentum)
 
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(200)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(20000)# (period)
dcd.unwrap(True)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(100000)# (period)
#xml.setOutputVelocity(True)
xml.setOutputMass(True)
xml.setOutputImage(True)
xml.setOutputBond(True)
xml.setOutputCharge(True)
app.add(xml)

#ready ro run
app.run(3000000)#(How many steps to run)
