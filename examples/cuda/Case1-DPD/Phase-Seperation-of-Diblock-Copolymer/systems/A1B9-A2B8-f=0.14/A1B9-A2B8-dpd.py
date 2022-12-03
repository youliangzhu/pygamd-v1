#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()
 
filename = 'A1B9-A2B8.xml' # initial configuration file
randomnum = 12340
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)# assign GPU by index
all_info = gala.AllInfo(build_method, perform_config) # build system information
 
dt = 0.04
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step

neighbor_list = gala.NeighborList(all_info, 1.0 ,0.05)#(,rcut,rbuffer)
dpd = gala.DPDForce(all_info, neighbor_list, 1.0, randomnum)#(,,rcut, the seed for RNG)
dpd.setParams('A', 'A', 25.0, 3.0)#(type,type,alpha,sigma)
dpd.setParams('A', 'B', 57.0, 3.0)#(type,type,alpha,sigma)
dpd.setParams('B', 'B', 25.0, 3.0)#(type,type,alpha,sigma)
app.add(dpd)
 
bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('A-A', 4.0, 0.0)#(bond type, K0-spring constant, R0-equilibrium distance)
bondforce.setParams('A-B', 4.0, 0.0)#(bond type, K0-spring constant, R0-equilibrium distance)
bondforce.setParams('B-B', 4.0, 0.0)#(bond type, K0-spring constant, R0-equilibrium distance)
app.add(bondforce)

group = gala.ParticleSet(all_info, "all") # a collection of particles
comp_info = gala.ComputeInfo(all_info, group)  # calculating system informations, such as temperature, pressure, and momentum

Gwvv = gala.DPDGWVV(all_info, group) # integration method with GWVV algorithm
app.add(Gwvv)

sort_method = gala.Sort(all_info)  # memory sorting to improve performance 
sort_method.setPeriod(80)
app.add(sort_method) 

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log') # output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(200)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')  # output the configuration files in mol2 format
mol2.setPeriod(0)# (period)
mol2.deleteBoundaryBond(True)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)  # output the configuration files in DCD format
dcd.setPeriod(10000)# (period)
dcd.unwrap(True)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')  # output the configuration files in xml format
xml.setPeriod(5000)# (period)
xml.setOutputBond(True)
app.add(xml)

#ready ro run
app.run(1600000)#(the number of steps to run)
neighbor_list.printStats() # output the information about neighbor_list 
