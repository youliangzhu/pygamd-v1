#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'Janus-one-patch-init.xml' # initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)  # assign GPU by index
all_info = gala.AllInfo(build_method,perform_config)# build system information

dt = 0.002
app = gala.Application(all_info, dt)  # build up an application with system information and integration time-step

# Build neighbor list object
neighbor_list = gala.NeighborList(all_info, 1.0 ,0.4) # rcut, rbuffer
neighbor_list.setBlockSize(96)
# neighbor_list.setDataReproducibility() #Data Reproducibility 

# build LZW force calcualte object
LZW = gala.LZWForce(all_info, neighbor_list, 1.0) # rcut
LZW.setParams('A', 'A' , 396.0, 1.0, 0.5, 242.0, 115.0/180.0*3.1415926)#(,,alpha,mu,nu,alpha1,beta)
LZW.setParams('A', 'B' , 396.0, 1.0, 0.5, 242.0, 115.0/180.0*3.1415926)#(,,alpha,mu,nu,alpha1,beta)
LZW.setParams('B', 'B' , 396.0, 1.0, 0.5, 242.0, 115.0/180.0*3.1415926)#(,,alpha,mu,nu,alpha1,beta)
LZW.setMethod('Janus')#Disk,Janus,ABAtriJanus,BABtriJanus
LZW.setBlockSize(96)
app.add(LZW)

group = gala.ParticleSet(all_info, "all")# a collection of particles
comp_info = gala.ComputeInfo(all_info, group)# calculating system informations, such as temperature, pressure, and momentum
Bere = gala.BerendsenAniNVT(all_info, group, comp_info, 1.0, 30.0*dt, 10.0*dt)#(,,temperature, tau, tauR)
app.add(Bere)
 
#v = gala.VariantLinear()
#v.setPoint(10000,1.0) # timesteps, temperature
#v.setPoint(20000,0.5)
#v.setPoint(30000,1.0)
#v.setPoint(40000,0.5)
#Bere.setT(v)

sort_method = gala.Sort(all_info) # sorting memory to improve performance 
sort_method.setPeriod(300)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log') # output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(10000)
DInfo.dumpAnisotropy()
app.add(DInfo)

# binary = gala.BinaryDump(all_info, 'particle')
# binary.setPeriod(100000)# (period)
# binary.setOutputForRestart()
# application.addStatistics(binary)

ZeroMomentum = gala.ZeroMomentum(all_info)  # remove the momentum of the center of mass
ZeroMomentum.setPeriod(10000)# (period)
app.add(ZeroMomentum)

xml = gala.XMLDump(all_info, 'particles') # output the configuration files in xml format
xml.setPeriod(500000)
xml.setOutputOrientation(True)
app.add(xml)

#for i in range(1,800)
#	application.run(10000)
#	xml.setPeriod(1)
#	application.run(10)
#	xml.setPeriod(10000)
app.run(16000000)# the number of time steps to run

neighbor_list.printStats() # output the information about neighbor_list 
#
#
