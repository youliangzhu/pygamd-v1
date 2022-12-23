#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'patchyinit.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig( int(_options.gpu))
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.002
app = gala.Application(all_info, dt)

# Build neighbor list object
neighbor_list = gala.NeighborList(all_info, 1.0 ,0.1)
neighbor_list.setBlockSize(96)
# neighbor_list.setDataReproducibility() #Data Reproducibility 

# build ani force calcualte object
ani = gala.AniForce(all_info, neighbor_list, 1.0)
ani.setParams('A', 'A' , 396.0, 2.0)#(,,alpha_R,mu)
ani.setParams('A', 'B' , 396.0, 2.0)#(,,alpha_R,mu)
ani.setParams('B', 'B' , 396.0, 2.0)#(,,alpha_R,mu)
ani.setPatches('patch-3.log')
ani.setBlockSize(96)
app.add(ani)


group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

bgroup = gala.ParticleSet(all_info, 'body')
rigidnvt = gala.NVTRigid(all_info, bgroup, 1.0, 0.2)
app.add(rigidnvt)

nbgroup = gala.ParticleSet(all_info,'non_body')
comp_info_nb = gala.ComputeInfo(all_info, nbgroup)
nh = gala.NoseHooverNVT(all_info, nbgroup, comp_info_nb, 1.0, 1.0)#( ,temperature, tau)
app.add(nh)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(300)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(10000)
DInfo.dumpAnisotropy()
app.add(DInfo)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(10000)# (period)
app.add(ZeroMomentum)


xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(500000)
xml.setOutputQuaternion(True)
xml.setOutputPatch(ani)
app.add(xml)

app.run(8000000)

neighbor_list.printStats()
#
#
