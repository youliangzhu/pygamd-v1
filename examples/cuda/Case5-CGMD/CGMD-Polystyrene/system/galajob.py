#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'coor.xml'# initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)# assign GPU by index
all_info = gala.AllInfo(build_method, perform_config)# build system information

dt = 0.00228
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step

neighbor_list = gala.NeighborList(all_info, 1.3, 0.15)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds() # exclude 1-2 interaction
neighbor_list.addExclusionsFromAngles()# exclude 1-3 interaction

#build pair force calculation object, 1 represent R monomer and 2 represent S monomer.
pair1 = gala.PairForceTable(all_info, neighbor_list,  2000) # (,,the number of data points)
pair1.setParams('1', '1' , 1.3 , "table.dat", 0, 3) #(type,type,rcut,potential data file,starting column,ending column)
pair1.setParams('1', '2' , 1.3 , "table.dat", 0, 3)#(type,type,rcut,potential data file,starting column,ending column)
pair1.setParams('2', '2' , 1.3 , "table.dat", 0, 3)#(type,type,rcut,potential data file,starting column,ending column)
pair1.setBlockSize(96)
app.add(pair1)

bond = gala.BondForceTable(all_info, 2000)# (,interaction length range,the number of data points)
bond.setParams('11', 2.0, "table.dat", 0, 3)#(bond type,potential data file,starting column,ending column)
bond.setParams('12', 2.0, "table.dat", 4, 7)#(bond type,potential data file,starting column,ending column)
bond.setBlockSize(96)
app.add(bond)

angle = gala.AngleForceTable(all_info, 500)# (,the number of data points)
angle.setParams('111' , "table.dat", 0, 3)#(angle type,potential data file,starting column,ending column)
angle.setParams('121' , "table.dat", 4, 7)#(angle type,potential data file,starting column,ending column)
angle.setParams('112' , "table.dat", 8, 11)#(angle type,potential data file,starting column,ending column)
angle.setParams('222' , "table.dat", 0, 3)#(angle type,potential data file,starting column,ending column)
angle.setParams('212' , "table.dat", 4, 7)#(angle type,potential data file,starting column,ending column)
angle.setParams('221' , "table.dat", 8, 11)#(angle type,potential data file,starting column,ending column)
angle.setParams('211' , "table.dat", 8, 11)#(angle type,potential data file,starting column,ending column)
angle.setParams('122' , "table.dat", 8, 11)#(angle type,potential data file,starting column,ending column)
angle.setBlockSize(96)
app.add(angle)

group = gala.ParticleSet(all_info, "all")# a collection of particles
comp_info = gala.ComputeInfo(all_info, group)  # calculating system informations, such as temperature, pressure, and momentum
Temperature = 500.0  #k
T = Temperature/300.0#reduced unit
P = 0.02446     #pressure 0.02446=101.325KPa
#NH = gala.NoseHoover(all_info, group, comp_info, T, 0.228)#( ,temperature, tau)
#app.add(NH)
npt = gala.NPT(all_info, group, comp_info, comp_info, T, P, 0.228, 2.28)# temperature,pressure,tauT,tauP
app.add(npt)
#Beren = gala.BerendsenNVT(all_info, group, comp_info, T, 0.228) 
#app.add(Beren) 

ZeroMomentum = gala.ZeroMomentum(all_info) # removing the momentum of the center of mass
ZeroMomentum.setPeriod(20000)# (period)
app.add(ZeroMomentum)

sort_method = gala.Sort(all_info) # memory sorting to improve data reading performance 
sort_method.setPeriod(300)# (period)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')# output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(1000)# (period)
app.add(DInfo)

#mol2 = gala.MOL2Dump(all_info, 'particles')
#mol2.setPeriod(0)# (period)
#app.add(mol2)

#dcd = gala.DCDDump(all_info, 'particle',True)
#dcd.setPeriod(10000)# (period)
#app.add(dcd)

xml = gala.XMLDump(all_info, 'particle')# output the configuration files in xml format
xml.setPeriod(10000)# (period)
xml.setOutputImage(True)
xml.setOutputType(True)
xml.setOutputBond(True)
xml.setOutputAngle(True)
xml.setOutputVelocity(True)
app.add(xml)

#binary = gala.BinaryDump(all_info, 'particle')
#binary.setPeriod(10000)# (period)
#binary.setOutputCtVersion(True)
#binary.setOutputImage(False)
#binary.setOutputVelocity(True)
#binary.setOutputType(False)
#binary.enableCompression(True)
#app.add(binary)

binary2 = gala.BinaryDump(all_info, 'particle')# output the configuration files in binary format
binary2.setPeriod(100000)# (period)
binary2.setOutputVelocity(True)
binary2.setOutputBond(True)
binary2.setOutputAngle(True)
binary2.setOutputForRestart()
app.add(binary2)

#ready ro run 
app.run (100000)#(the number of steps to run)
neighbor_list.printStats()# output the information about neighbor_list 
#
#
