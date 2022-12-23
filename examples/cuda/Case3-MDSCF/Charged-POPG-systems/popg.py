#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options
 
filename = 'popg.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)
 
epsilonr=80.0
dt = 0.03
Temperature = 301.15  #k
T = Temperature*8.3143/1000.0#reduced unit
app = gala.Application(all_info, dt)

scf = gala.MDSCFForce(all_info,    32,    32,    36, 0.100)#(mx,my,mz,compressibility)
scf.setParams('L', 'L',  0.0)#(chi)
scf.setParams('L', 'P', -3.6)#(chi)
scf.setParams('L', 'G',  4.5)#(chi)
scf.setParams('L', 'C',  13.25)#(chi)
scf.setParams('L', 'D',  9.3)#(chi)
scf.setParams('L', 'Na', -6.0)#(chi)
scf.setParams('L', 'W',  0.0)#(chi)

scf.setParams('P', 'P',  0.0)#(chi)
scf.setParams('P', 'G',  2.5)#(chi)
scf.setParams('P', 'C', 13.5)#(chi)
scf.setParams('P', 'D', 11.7)#(chi)
scf.setParams('P', 'Na', 0.0)#(chi)
scf.setParams('P', 'W', -3.6)#(chi)

scf.setParams('G', 'G',  0.0)#(chi)
scf.setParams('G', 'C',  8.3)#(chi)
scf.setParams('G', 'D',  8.3)#(chi)
scf.setParams('G', 'Na', -3.0)#(chi)
scf.setParams('G', 'W',  4.5)#(chi)

scf.setParams('C', 'C',  0.0)#(chi)
scf.setParams('C', 'D',  0.0)#(chi)
scf.setParams('C', 'Na', 13.5)#(chi)
scf.setParams('C', 'W',  33.75)#(chi)

scf.setParams('D', 'D',  0.0)#(chi)
scf.setParams('D', 'Na', 9.3)#(chi)
scf.setParams('D', 'W', 23.25)#(chi)

scf.setParams('Na', 'Na', 0.0)#(chi)
scf.setParams('Na', 'W',  0.0)#(chi)

scf.setParams('W', 'W',   0.0)#(chi)
scf.setPeriodScf(1,  100)#(idl2_step, idl_step)
app.add(scf)

pfme = gala.PFMEForce(all_info,    32,   32,   36,  2.45,  epsilonr)#(mx,my,mz, kappa, epsilonr)
pfme.setPeriodPFME(1,  100)#(idl2_step, idl_step)
app.add(pfme)

bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('L-P',  1250.000,  0.370)#(,K0, R0)
bondforce.setParams('P-G',  1250.000,  0.370)#(,K0, R0)
bondforce.setParams('G-G',  1250.000,  0.370)#(,K0, R0)
bondforce.setParams('G-C',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('C-C',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('C-D',  1250.000,  0.470)#(,K0, R0)
app.add(bondforce)
 
angleforce = gala.AngleForceHarmonicCos(all_info)
angleforce.setParams('P-G-G',   25.000, 120.000)#(,K0, t0(degree))
angleforce.setParams('P-G-C',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('G-C-C',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('C-C-C',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('G-C-D',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('C-D-C',   45.000, 120.000)#(,K0, t0(degree))
angleforce.setParams('C-C-D',   25.000, 180.000)#(,K0, t0(degree))
app.add(angleforce)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

thermo = gala.AndersenNVT(all_info, group, T, 7.000, 1234)#( ,temperature, coll_freq, random seed)
app.add(thermo)

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
xml.setPeriod(1000000)# (period)
xml.setOutputBond(True)
xml.setOutputAngle(True)
xml.setOutputImage(True)
xml.setOutputCharge(True)
app.add(xml)

#ready ro run
app.run( 20000000)#(How many steps to run)

