#!/usr/bin/python2
from poetry import cu_gala as gala 
from poetry import _options

filename = 'slj.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.001
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 3.0, 0.1)#(,rcut,rbuffer)
neighbor_list.setFilterDiameters()

lj = gala.SLJForce(all_info, neighbor_list, 3.0)
lj.setParams('A', 'A' ,1.0 ,1.0 ,1.0)
lj.setParams('A', 'B' ,1.0 ,1.0 ,1.0)
lj.setParams('B', 'B' ,1.0 ,1.0 ,1.0)
lj.setEnergy_shift()
app.add(lj)

group = gala.ParticleSet(all_info,'all')
comp_info = gala.ComputeInfo(all_info, group)
nvt = gala.NoseHooverNVT(all_info, group, comp_info, 1.0, 1.0)
app.add(nvt)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(100)
app.add(DInfo)

zm = gala.ZeroMomentum(all_info)
zm.setPeriod(100000)
app.add(zm)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(200)# (period)
app.add(sort_method)

mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)

dcd = gala.DCDDump(all_info, 'particle',True)
dcd.setPeriod(1000)# (period)
dcd.unwrap(True)
app.add(dcd)

xml = gala.XMLDump(all_info, 'particle')
xml.setPeriod(10000)# (period)
xml.setOutputImage(True)
xml.setOutputVelocity(True)
xml.setOutputDiameter(True)
xml.setOutputType(True)
app.add(xml)

#ready ro run 
app.run(20000)
neighbor_list.printStats()
#
#

