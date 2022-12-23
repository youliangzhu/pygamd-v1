#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'example.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.001
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 3.0 ,0.4)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()
neighbor_list.addExclusionsFromAngles()
neighbor_list.addExclusionsFromDihedrals()

lj = gala.LJForce(all_info, neighbor_list, 3.0)#(,,rcut)
lj.setParams('A', 'A' ,1.0 ,1.0 ,1.0)#(,,epsilon,sigma,alpha)
app.add(lj)

bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('A-A',  1000.0,  0.5)#(,K0, R0)
app.add(bondforce)
 
angleforce = gala.AngleForceHarmonic(all_info)
angleforce.setParams('A-A-A', 100.0, 120.0)#(,K0, t0(degree))
app.add(angleforce)

dihedralforce = gala.DihedralForceHarmonic(all_info)
dihedralforce.setParams('A-A-A-A', 25.0, 180.0)#(,K0, t0(degree))
app.add(dihedralforce)


group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

nh = gala.NoseHooverNVT(all_info, group, comp_info, 1.0, 0.5)#( ,temperature, tau)
app.add(nh)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(1000)# (period)
app.add(sort_method)

xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(10000)# (period)
app.add(xml)


dInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
dInfo.setPeriod(300)# (period)
app.add(dInfo)

#ready ro run 
app.run(20000)
neighbor_list.printStats()
#
#
