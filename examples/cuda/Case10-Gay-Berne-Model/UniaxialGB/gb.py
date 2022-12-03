#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()

filename = 'gb.xml' # initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)# assign GPU by index
all_info = gala.AllInfo(build_method, perform_config)# build system information

dt = 0.001
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step

# Build neighbor list object
neighbor_list = gala.NeighborList(all_info, 10.0 ,1.0)#(,rcut,rbuffer)
neighbor_list.setBlockSize(96)
# neighbor_list.setDataReproducibility() #Data Reproducibility 

nu=1.0
mu=2.0
sigmae=3.0
sigmas=1.0
epsilone=0.5
epsilons=3.0
ps=1.0
# build GB force calculation object
GB = gala.GBForce(all_info, neighbor_list, 10.0)#(,,rcut)
GB.setParams('A', 'A' , 1.5, 1.5, nu, mu, sigmae, sigmas, epsilone, epsilons, ps)#(type,type,epsilon0,sigma0,nu,mu,sigmae,sigmas,epsilone,epsilons, ps)
GB.setBlockSize(96)
app.add(GB)

group = gala.ParticleSet(all_info, "all")# a collection of particles
comp_info = gala.ComputeInfo(all_info, group)# calculating system informations, such as temperature, pressure, and momentum
Bere = gala.BerendsenAniNVT(all_info, group, comp_info, 1.0, 30.0*dt, 10.0*dt)#(,,temperature, tau, tauR)
app.add(Bere)

sort_method = gala.Sort(all_info) # sorting memory to improve performance 
sort_method.setPeriod(300)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')  # memory sorting to improve performance
DInfo.setPeriod(1000)
DInfo.dumpAnisotropy()
app.add(DInfo)

xml = gala.XMLDump(all_info, 'particles')# output the configuration files in xml format
xml.setPeriod(10000)
xml.setOutputOrientation(True)
xml.setOutputEllipsoid(GB)
app.add(xml)

app.run(1000000)#(the number of steps to run)
neighbor_list.printStats()# output the information about neighbor_list 
#
#
