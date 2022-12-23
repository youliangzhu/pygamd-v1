#!/usr/bin/python
from poetry import force_field_gala
from poetry import cu_gala as gala 
from poetry import _options
 
filename = 'Equ.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)
 
dt = 0.03
app = gala.Application(all_info, dt)

const = 138.935
epsilon_r = 15

neighbor_list = gala.NeighborList(all_info, 1.2, 0.12)#(,rcut,rbuffer)
neighbor_list.exclusion(["bond", "constraint", "dihedral", "vsite"])

#intramolecular exclusion for cholesterol
offset=25272
nmol=900
for i in range(0,nmol):
	idx=offset+i*8
	for j in range(0, 8):
		for k in range(j+1, 8):		
			neighbor_list.addExclusion(idx+j,idx+k)

lj = force_field_gala.LJCoulombShiftForce(all_info, neighbor_list, 1.2, 0.9, epsilon_r, "Equ.force_field")
app.add(lj)

bondforce = force_field_gala.BondForceHarmonic(all_info, "Equ.force_field")
app.add(bondforce)

angleforce = force_field_gala.AngleForceHarmonicCos(all_info, "Equ.force_field")
app.add(angleforce)

dihedralforce = force_field_gala.DihedralForceHarmonic(all_info, "Equ.force_field")
app.add(dihedralforce)

bond_constraint = force_field_gala.BondConstraint(all_info, "Equ.force_field")
bond_constraint.setExpansionOrder(4)
bond_constraint.setNumIters(1)
app.add(bond_constraint)

vs = force_field_gala.Vsite(all_info, "Equ.force_field")
app.add(vs)

group = gala.ParticleSet(all_info, "all")
comp_info = gala.ComputeInfo(all_info, group)

Temperature = 310.000  #k
T = Temperature*8.3143/1000.0#reduced unit
thermo = gala.NPTMTKSD(all_info, group, comp_info, comp_info, T, 0.0, 1.0, 12.0)
thermo.setSemiisotropic(0.0, 0.0)
thermo.setCompressibility(3e-4, 3e-4, 0.0)
app.add(thermo)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(800)
app.add(sort_method)

#ZeroMomentum = gala.ZeroMomentum(all_info)
#ZeroMomentum.setPeriod(100000)# (period)
#app.add(ZeroMomentum)
 
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(400)
DInfo.dumpVirialEnergy(lj)
app.add(DInfo)

mol2 = gala.MOL2Dump(all_info, 'particles')
mol2.setPeriod(0)# (period)
app.add(mol2)
 
dcd = gala.DCDDump(all_info, 'particles',True)
dcd.setPeriod(50000)# (period)
app.add(dcd)
 
xml = gala.XMLDump(all_info, 'particles')
xml.setPeriod(10000)# (period)
xml.setOutput(["velocity", "image", "charge",  "mass", "bond", "constraint", "vsite", "angle", "dihedral"])
app.add(xml)

#ready to run
app.run(1)#(How many steps to run)

