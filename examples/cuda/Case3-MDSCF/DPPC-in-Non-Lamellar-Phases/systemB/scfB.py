#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()
 
filename = 'scfB.xml'# initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)# assign GPU by index
all_info = gala.AllInfo(build_method,perform_config)# build system information
 
dt = 0.03
Temperature = 318.000  #k
T = Temperature*8.3143/1000.0#reduced unit
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step

# intra-molecular neighbor list
neighbor_list = gala.IntraMolList(all_info, 0.9000, 0.3000)#(rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()

# intra-molecular interactions
lj = gala.LJForce(all_info, neighbor_list, 0.9000)
lj.setParams('N', 'N',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('N', 'P',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('N', 'G',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('N', 'C',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('N', 'W',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('P', 'P',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('P', 'G',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('P', 'C',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('P', 'W',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('G', 'G',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('G', 'C',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('G', 'W',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('C', 'C',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('C', 'W',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
lj.setParams('W', 'W',  2.000,  0.470,  1.000)#(type,type,epsilon,sigma,alpha)
app.add(lj)

# MD-SCF inter-molecular force
scf = gala.MDSCFForce(all_info,    22,    22,    22, 0.100)#(mx,my,mz,compressibility)
scf.setParams('N', 'N',  0.000)#(chi)
scf.setParams('N', 'P', -1.500)#(chi)
scf.setParams('N', 'G',  6.300)#(chi)
scf.setParams('N', 'C',  9.000)#(chi)
scf.setParams('N', 'W', -8.130)#(chi)
scf.setParams('P', 'P',  0.000)#(chi)
scf.setParams('P', 'G',  4.500)#(chi)
scf.setParams('P', 'C', 13.500)#(chi)
scf.setParams('P', 'W', -3.600)#(chi)
scf.setParams('G', 'G',  0.000)#(chi)
scf.setParams('G', 'C',  6.300)#(chi)
scf.setParams('G', 'W',  4.500)#(chi)
scf.setParams('C', 'C',  0.000)#(chi)
scf.setParams('C', 'W', 33.750)#(chi)
scf.setParams('W', 'W',  0.000)#(chi)
scf.setPeriodScf(   1,  300)#(idl2_step, idl_step)
app.add(scf)

# bond stretching interaction by harmonic potential
bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('N-P',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('P-G',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('G-G',  1250.000,  0.370)#(,K0, R0)
bondforce.setParams('G-C',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('C-C',  1250.000,  0.470)#(,K0, R0)
app.add(bondforce)

# angle bending interaction 
angleforce = gala.AngleForceHarmonicCos(all_info)
angleforce.setParams('P-G-G',   25.000, 120.000)#(,K0, t0(degree))
angleforce.setParams('P-G-C',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('G-C-C',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('C-C-C',   25.000, 180.000)#(,K0, t0(degree))
app.add(angleforce)

group = gala.ParticleSet(all_info, "all")# a collection of particles
comp_info = gala.ComputeInfo(all_info, group)  # calculating system informations, such as temperature, pressure, and momentum
thermo = gala.AndersenNVT(all_info, group, T, 7.000, 1234)#( ,temperature, coll_freq, random seed)
app.add(thermo)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')# output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(200)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')# output the configuration files in mol2 format
mol2.setPeriod(100000)# (period)
mol2.deleteBoundaryBond(True)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)# output the configuration files in DCD format
dcd.setPeriod(10000)# (period)
dcd.unwrap(True)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles') # output the configuration files in xml format
xml.setPeriod(100000)# (period)
xml.setOutputBond(True)
xml.setOutputImage(True)
app.add(xml)

#ready ro run
app.run( 20000000)#(the number of steps to run)

