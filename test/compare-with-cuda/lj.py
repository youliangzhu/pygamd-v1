#!/usr/bin/python
from poetry import cu_gala as gala
from poetry import _options
 

filename = 'lj.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.001
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 3.0 ,0.4)#(,rcut,rbuffer)
lj = gala.LJForce(all_info, neighbor_list, 3.0)#(,,rcut)
lj.setParams('a', 'a' ,1.0 ,1.0 ,1.0)#(,,epsilon,sigma,alpha)
app.add(lj)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

Nh = gala.NoseHooverNVT(all_info, group, comp_info, 1.0, 0.5)#( ,temperature, tau)
app.add(Nh)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(1000)# (period)
app.add(sort_method)

xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(10000)# (period)
app.add(xml)


DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(300)# (period)
app.add(DInfo)

#ready ro run 
app.run(2000000)
neighbor_list.printStats()
#
#
