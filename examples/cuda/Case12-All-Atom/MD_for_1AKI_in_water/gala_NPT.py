#!/usr/bin/python
from poetry import force_field_itp
from poetry import cu_gala as gala 
from poetry import _options

filename = 'eq_NPT.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.001
app = gala.Application(all_info, dt)
neighbor_list = gala.NeighborList(all_info, 1.0, 0.1)#(,rcut,rbuffer)
neighbor_list.exclusion(["bond", "angle", "dihedral"])

lj = force_field_itp.LJEwaldForce(all_info, neighbor_list, 1.0, "ffnonbonded.itp")
lj.setEnergy_shift()
app.add(lj)

groupC = gala.ParticleSet(all_info, "charge")
pppm = gala.PPPMForce(all_info, neighbor_list, groupC)
pppm.setParams(0.15, 4, 1.0) #grid space, spread order, rcut in real space
pppm.setChargeCellList(False)
app.add(pppm)

bondforce = force_field_itp.BondForceHarmonic(all_info, "ffbonded.itp")
app.add(bondforce)

angleforce = force_field_itp.AngleForceHarmonic(all_info, "ffbonded.itp")
app.add(angleforce)

dihedralforce = force_field_itp.DihedralForceAmberCosine(all_info, "ffbonded.itp")
app.add(dihedralforce)

pairforce = gala.LJCoulombPair(all_info)
pairforce.setDividedFactorVDWELEC('p1', 2, 1.200048002)
app.add(pairforce)

Temperature = 298.000  #k
T = Temperature*8.3143/1000.0#reduced unit

Pressure = 1.0/16.3882449645417 #bar
P = Pressure

group = gala.ParticleSet(all_info,'all')
comp_info = gala.ComputeInfo(all_info, group)
npt = gala.NPTMTK(all_info, group, comp_info, comp_info, T, P, 0.2, 2.0)
npt.setSemiisotropic(0.1, 0.1)
app.add(npt)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(1000)# (period)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(1000)# (period)
app.add(ZeroMomentum)


DInfo = gala.DumpInfo(all_info, comp_info, 'data_NPT.log')
DInfo.setPeriod(1000)# (period)
DInfo.dumpPotential(lj)
DInfo.dumpPotential(pppm)
DInfo.dumpPotential(bondforce)
DInfo.dumpPotential(angleforce)
DInfo.dumpPotential(dihedralforce)
DInfo.dumpPotential(pairforce)
DInfo.dumpBoxSize()
app.add(DInfo)


xml = gala.XMLDump(all_info, 'NPT')
xml.setOutput(['image', 'bond'])
xml.setPeriod(100000)
app.add(xml)


# Period1
#ready ro run 
app.run(200000)
neighbor_list.printStats()

