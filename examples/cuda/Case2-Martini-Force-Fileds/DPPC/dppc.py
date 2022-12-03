#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()
 
filename = 'dppc.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)
 
dt = 0.02
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 1.2, 0.12)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()
epsilon_r=15.0

pf = gala.LJCoulombShiftForce(all_info, neighbor_list)
pf.setParams('P4', 'P4', 4.99986405887, 0.470000564947, 1.0, 1.2, 0.9)
pf.setParams('Na', 'Na', 3.99979176977, 0.470000696296, 1.0, 1.2, 0.9)
pf.setParams('C1', 'C1', 3.50000431856, 0.47000041484,  1.0, 1.2, 0.9)
pf.setParams('Qa', 'Qa', 4.99986405887, 0.470000564947, 1.0, 1.2, 0.9)
pf.setParams('Q0', 'Q0', 3.50000431856, 0.47000041484,  1.0, 1.2, 0.9)
pf.setParams('P4', 'Na', 3.99979176977, 0.470000696296, 1.0, 1.2, 0.9)
pf.setParams('P4', 'C1', 1.99997049288, 0.470000499274, 1.0, 1.2, 0.9)
pf.setParams('P4', 'Qa', 5.59977163138, 0.470001759455, 1.0, 1.2, 0.9)
pf.setParams('P4', 'Q0', 5.59977163138, 0.470001759455, 1.0, 1.2, 0.9)
pf.setParams('Na', 'C1', 2.70013873615, 0.469998985738, 1.0, 1.2, 0.9)
pf.setParams('Na', 'Qa', 3.99979176977, 0.470000696296, 1.0, 1.2, 0.9)
pf.setParams('Na', 'Q0', 3.99979176977, 0.470000696296, 1.0, 1.2, 0.9)
pf.setParams('C1', 'Qa', 1.99999380085, 0.619999891705, 1.0, 1.2, 0.9)
pf.setParams('C1', 'Q0', 1.99999380085, 0.619999891705, 1.0, 1.2, 0.9)
pf.setParams('Qa', 'Q0', 4.49982791432, 0.470000623324, 1.0, 1.2, 0.9)
pf.setCoulomb(1.2, 0.9, epsilon_r)
app.add(pf)

bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('Q0-Qa',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('Qa-Na',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('Na-Na',  1250.000,  0.370)#(,K0, R0)
bondforce.setParams('Na-C1',  1250.000,  0.470)#(,K0, R0)
bondforce.setParams('C1-C1',  1250.000,  0.470)#(,K0, R0)
app.add(bondforce)
 
angleforce = gala.AngleForceHarmonicCos(all_info)
angleforce.setParams('Qa-Na-Na',   25.000, 120.000)#(,K0, t0(degree))
angleforce.setParams('Qa-Na-C1',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('Na-C1-C1',   25.000, 180.000)#(,K0, t0(degree))
angleforce.setParams('C1-C1-C1',   25.000, 180.000)#(,K0, t0(degree))
app.add(angleforce)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

Temperature = 318.000  #k
T = Temperature*8.3143/1000.0#reduced unit
thermo = gala.NoseHooverNVT(all_info, group, comp_info, T, 0.5)
app.add(thermo)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(400)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(100000)# (period)
app.add(ZeroMomentum)
 
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(400)
app.add(DInfo)
 
mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(100000)# (period)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(100000)# (period)
xml.setOutputBond(True)
xml.setOutputAngle(True)
xml.setOutputVelocity(True)
xml.setOutputMass(True)
app.add(xml)

#ready ro run

app.run(  1000000)#(How many steps to run)

