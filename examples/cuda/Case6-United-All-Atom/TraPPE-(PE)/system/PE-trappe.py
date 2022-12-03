#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()
 
filename = 'PE-1000.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)
 
dt = 0.002
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 1.4, 0.1)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()
neighbor_list.addExclusionsFromAngles()
neighbor_list.addExclusionsFromDihedrals()

pf = gala.LJForce(all_info, neighbor_list, 1.4)
pf.setParams('C_33', 'C_33',  0.814818,  0.375, 1.0)#(,,epsilon,sigma,alpha,rcut)
pf.setParams('C_33', 'C_32',  0.558246,  0.385, 1.0)#(,,epsilon,sigma,alpha,rcut)
pf.setParams('C_32', 'C_32',  0.382465,  0.395, 1.0)#(,,epsilon,sigma,alpha,rcut)
app.add(pf)


bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('C_33-C_32',  129704,  0.154)#(,K0, R0)
bondforce.setParams('C_32-C_32',  129704,  0.154)#(,K0, R0)
app.add(bondforce)
 
angleforce = gala.AngleForceHarmonic(all_info)
angleforce.setParams('C_33-C_32-C_32', 519.6545, 114.0)#(,K0, t0(degree))
angleforce.setParams('C_32-C_32-C_32', 519.6545, 114.0)#(,K0, t0(degree))
app.add(angleforce)

dihedralforce = gala.DihedralForceOPLSCosine(all_info)
dihedralforce.setParams('C_33-C_32-C_32-C_32', 0.0, 2.95188, -0.566963, 6.57940, 0.0)#(,K0, t0(degree))
dihedralforce.setParams('C_32-C_32-C_32-C_32', 0.0, 2.95188, -0.566963, 6.57940, 0.0)#(,K0, t0(degree))
app.add(dihedralforce)


group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)


Temperature0 = 450.000  #k
Temperature1 = 300.000  #k
T0 = Temperature0*8.3143/1000.0#reduced unit
T1 = Temperature1*8.3143/1000.0#reduced unit

T=gala.VariantLinear()
T.setPoint(0,T0)
T.setPoint(10000000,T1)
thermo = gala.NPT(all_info, group, comp_info,comp_info, T0, 0.06102, 1.0, 1.0)
thermo.setT(T)
app.add(thermo)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(400)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(0)# (period)
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
xml.setOutputDihedral(True)
xml.setOutputVelocity(True)
xml.setOutputMass(True)
xml.setOutputImage(True)
app.add(xml)

#ready ro run
app.run(  100000000)#(How many steps to run)

