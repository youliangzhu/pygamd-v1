#!/usr/bin/python
from poetry import force_field_itp
from poetry import cu_gala as gala 
from poetry import _options

filename = 'conf.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.002
app = gala.Application(all_info, dt)
neighbor_list = gala.NeighborList(all_info, 1.0, 0.1)#(,rcut,rbuffer)
neighbor_list.exclusion(["bond", "angle", "dihedral"])

lj = force_field_itp.LJEwaldForce(all_info, neighbor_list, 1.0, "ffnonbonded.itp")
lj.setEnergy_shift()
lj.setDispVirialCorr(True)#dispersion virial correction
app.add(lj)

groupC = gala.ParticleSet(all_info, "charge")
pppm = gala.PPPMForce(all_info, neighbor_list, groupC)
pppm.setParams(0.12, 4, 1.0) #grid space, spread order, rcut in real space
pppm.setChargeCellList(False)
app.add(pppm)

bondforce = force_field_itp.BondForceHarmonic(all_info, "ffbonded.itp")
app.add(bondforce)

angleforce = force_field_itp.AngleForceHarmonic(all_info, "ffbonded.itp")
app.add(angleforce)

dihedralforce = force_field_itp.DihedralForceAmberCosine(all_info, "ffbonded.itp")
app.add(dihedralforce)

Temperature = 298.000  #k
T = Temperature*8.3143/1000.0#reduced unit

group= gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)
npt = gala.NPTMTK(all_info, group, comp_info, comp_info, T, 0.0602, 1.0, 1.0)
#npt.setCompressibility(4.5e-5,4.5e-5,4.5e-5)
app.add(npt)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(800)# (period)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(100000)# (period)
app.add(ZeroMomentum)


DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(1000)# (period)
DInfo.dumpPotential(lj)
DInfo.dumpPotential(pppm)
DInfo.dumpPotential(bondforce)
DInfo.dumpPotential(angleforce)
DInfo.dumpPotential(dihedralforce)
app.add(DInfo)

xml = gala.XMLDump(all_info, 'p')
xml.setPeriod(10000)# (period)
xml.setOutput(["bond", "angle", "dihedral", "mass", "velocity", "charge","image"])
app.add(xml)

# Period1
#ready ro run 
app.run(1000000)
neighbor_list.printStats()

