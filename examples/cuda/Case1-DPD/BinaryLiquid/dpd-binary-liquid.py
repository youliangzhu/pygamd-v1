#!/usr/bin/python
#### script head #####
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()
 
#### read date from file #####
filename = 'AB.xml'
randomnum = 12340
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)
 
#### build application #####
dt = 0.01
T = 1.0             #temperature reduced unit
app = gala.Application(all_info, dt)


#### interactions #####
neighbor_list = gala.NeighborList(all_info, 1.0 ,0.05)#(,rcut,rbuffer)
dpd = gala.DPDForce(all_info, neighbor_list, 1.0, randomnum)#(,,rcut,seed)
dpd.setParams('A', 'A', 25.0, 3.0)#(,,p0,p1,rcut,func)
dpd.setParams('A', 'B', 40.0, 3.0)#(,,p0,p1,rcut,func)
dpd.setParams('B', 'B', 25.0, 3.0)#(,,p0,p1,rcut,func)
app.add(dpd) 

#### integration #####
group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

gwvv = gala.DPDGWVV(all_info, group)
app.add(gwvv)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(80)
app.add(sort_method) 

#### data dump #####
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(200)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(100000)# (period)
mol2.deleteBoundaryBond(True)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(10000)# (period)
dcd.unpbc(True)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(5000)# (period)
app.add(xml)

#### ready ro run ####
app.run( 50000)#(How many steps to run)
neighbor_list.printStats()
