#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'relax.xml'# initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)  # assign GPU by index
all_info = gala.AllInfo(build_method,perform_config)# build system information

dt = 0.001
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step

neighbor_list = gala.NeighborList(all_info, 3.0 ,0.4)#(,rcut,rbuffer)
lj = gala.LJForce(all_info, neighbor_list, 3.0)#(,,rcut)
lj.setParams("A", "A", 1.0, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
lj.setParams("B", "B", 1.0, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
lj.setParams("C", "C", 1.0, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
lj.setParams("A", "B", 0.5, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
lj.setParams("A", "C", 1.0, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
lj.setParams("B", "C", 1.0, 1.0, 1.0)#(type,type,epsilon,sigma,alpha)
app.add(lj)

bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('C-A', 40.0 , 1.0)#(,K0, R0)
bondforce.setParams('A-A', 40.0 , 1.0)#(,K0, R0)
bondforce.setParams('C-B', 40.0 , 1.0)#(,K0, R0)
bondforce.setParams('B-B', 40.0 , 1.0)#(,K0, R0)
app.add(bondforce)

group = gala.ParticleSet(all_info, 'all')
comp_info = gala.ComputeInfo(all_info, group)

T=1.8
Bd = gala.LangevinNVT(all_info, group, T, 123)#(,,temperature, seed)
Bd.setGamma('A', 1.0)#(,gama)
Bd.setGamma('B', 1.0)#(,gama)
Bd.setGamma('C', 1.0)#(,gama)
app.add(Bd)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(300)
app.add(sort_method) 

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(1000)
app.add(DInfo)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(5000)# (period)
app.add(ZeroMomentum)

mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(200000)# (period)
dcd.unwrap(True)
app.add(dcd)

#write bin file
binary2 = gala.BinaryDump(all_info, 'particle')
binary2.setPeriod(10000)# (period)
binary2.setOutputForRestart()
app.add(binary2)

xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(100000)# (period)
xml.setOutputBond(True)
xml.setOutputImage(True)
app.add(xml)

#ready ro run
app.run(20000)#(How many steps to run)
app.setDt(0.01)
app.run(20000000)#(How many steps to run)
neighbor_list.printStats()

