#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'particles.xml' # initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu) # GPU index
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.005
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 3.0 ,0.1)#(,rcut, rbuffer)
neighbor_list.addExclusionsFromBodys() # remove the interactions of particles in a same body

lj = gala.LJForce(all_info, neighbor_list,3.0)
lj.setParams('A', 'A' ,1.0 ,1.0 ,1.0, 1.12246) #type, type, epsilon, sigma, alpha, rcut
lj.setParams('A', 'B' ,1.0 ,1.0 ,1.0, 1.12246)#type, type, epsilon, sigma, alpha, rcut
lj.setParams('B', 'B' ,1.0 ,1.0 ,1.0, 3)#type, type, epsilon, sigma, alpha, rcut
lj.setEnergy_shift()#The potential is shifted to make the energy be zero at cutoff distance
app.add(lj)

# Brownian dynamics for rigid bodies
bgroup = gala.ParticleSet(all_info, 'body') # the group of the body particles
bdnvt_rigid = gala.BDNVTRigid(all_info, bgroup, 0.5, 123) # temperature, seed for random number generator 
app.add(bdnvt_rigid)

group = gala.ParticleSet(all_info,'all') # all particles 
comp_info = gala.ComputeInfo(all_info, group) # calculating system informations, such as temperature, pressure, and momentum

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log') # output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(100)
app.add(DInfo)

zm = gala.ZeroMomentum(all_info) # remove the momentum of the center of mass
zm.setPeriod(100) # period
app.add(zm)

sort_method = gala.Sort(all_info) # sorting memory to improve performance 
sort_method.setPeriod(200)# (period)
app.add(sort_method)

mol2 = gala.MOL2Dump(all_info, 'particles') # output the configuration files in mol2 format
mol2.setPeriod(0)# (period)
app.add(mol2)

dcd = gala.DCDDump(all_info, 'particle',True)# output the configuration files in DCD format
dcd.setPeriod(10000000)# (period)
#dcd.unwrap(True)
app.add(dcd)

xml = gala.XMLDump(all_info, 'particle') # output the configuration files in xml format
xml.setPeriod(10000000)# (period)
xml.setOutputImage(True)
xml.setOutputBond(True)
xml.setOutputVelocity(True)
xml.setOutputDiameter(True)
xml.setOutputType(True)
app.add(xml)

#ready ro run 
app.run(100000000) # the number of time steps to run
neighbor_list.printStats() # output the information about neighbor_list 
#
#
