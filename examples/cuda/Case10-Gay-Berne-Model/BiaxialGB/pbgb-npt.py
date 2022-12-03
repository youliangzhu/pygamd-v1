#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()

filename = 'gbinit.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.002
app = gala.Application(all_info, dt)

# Build neighbor list object
neighbor_list = gala.NeighborList(all_info, 4.0, 0.8)
# neighbor_list.setDataReproducibility() #Data Reproducibility 
neighbor_list.addExclusionsFromBodys()

# build ani force calcualte object
pbgb = gala.PBGBForce(all_info, neighbor_list)
pbgb.setGUM(1.0, 3.0, 1.0);#(gamma, niu, miu)
pbgb.setParams('B', 'B' , 2.0, 1.0, 4.0)#(,,epsilon, sigma, rcut)
pbgb.setParams('A', 'A' , 2.0, 1.0, 2.5)#(,,epsilon, sigma, rcut)
pbgb.setParams('A', 'B' , 1.0, 1.0, 3.25)#(,,epsilon, sigma, rcut)
pbgb.setAspheres('patch-2.log')#(,a,b,c,eia_one,eib_one,eic_one)
#pbgb.setPatches('patch-2.log')
app.add(pbgb)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

# bgroup = gala.ParticleSet(all_info, 'body')
# rigidnvt = gala.NVTRigid(all_info, bgroup, 0.6, 1.0)
# app.add(rigidnvt)

# nbgroup = gala.ParticleSet(all_info,'non_body')
# comp_info_nb = gala.ComputeInfo(all_info, nbgroup)
# Nh = gala.NoseHooverNVT(all_info, nbgroup, comp_info_nb, 0.6, 1.0)#( ,temperature, tau)
# app.add(Nh)


# build rigid nve integration
bgroup = gala.ParticleSet(all_info, 'body')
comp_info_b = gala.ComputeInfo(all_info, bgroup)
rigidnpt = gala.NPTRigid(all_info, bgroup, comp_info_b, comp_info, 2.4, 8.0, 0.5, 0.5)  # T P tau tauP
app.add(rigidnpt)


nbgroup = gala.ParticleSet(all_info,'non_body')
comp_info_nb = gala.ComputeInfo(all_info, nbgroup)
npt = gala.NPT(all_info, nbgroup, comp_info_nb, comp_info, 2.4, 8.0, 0.5, 0.5)# T P tau tauP
app.add(npt)


sort_method = gala.Sort(all_info)
sort_method.setPeriod(300)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(500)
DInfo.dumpAnisotropy()
app.add(DInfo)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(50000)# (period)
app.add(ZeroMomentum)


xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(20000)
xml.setOutputQuaternion(True)
xml.setOutputBody(True)
xml.setOutputMass(True)
xml.setOutputImage(True)
xml.setOutputEllipsoid(pbgb)
app.add(xml)

app.run(10000000)
neighbor_list.printStats()
#
#
