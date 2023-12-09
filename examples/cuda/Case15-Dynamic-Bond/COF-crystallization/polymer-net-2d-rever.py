#!/usr/bin/python
from poetry import cu_gala as gala
from poetry import _options


filename = 'pn2d.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.0005
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 2**(1.0/6.0), 0.1)#(,rcut,rbuffer)
#neighbor_list.addExclusionsFromBonds()
#neighbor_list.addExclusionsFromAngles()

lj = gala.LJForce(all_info, neighbor_list, 2**(1.0/6.0))# force rcut
lj.setParams("A", "A", 1.0, 1.0, 1.0) # type1, type2, epsilon, sigma, alpha, rcut
lj.setParams('B', 'B', 1.0, 1.0, 1.0)
lj.setParams('A', 'B', 1.0, 1.0, 1.0)
app.add(lj)

all_info.addBondTypeByPairs()
all_info.addAngleTypeByPairs()

bf = gala.BondForceHarmonic(all_info)
bf.setParams('A-B', 1250.0, 1.0)#(,K0, R0)
bf.setParams('B-B', 1250.0, 1.0)#(,K0, R0)
app.add(bf)

af = gala.AngleForceHarmonic(all_info)
af.setParams('B-A-B', 1000.0, 120)#(,K0, R0)
af.setParams('A-B-B', 1000.0, 180)#(,K0, R0)
app.add(af)

group = gala.ParticleSet(all_info,'all')
comp_info = gala.ComputeInfo(all_info, group)

T = 1.0                                    #reduced unit
bd=gala.LangevinNVT(all_info, group, T, 123) # all_info, group, T, seed
app.add(bd)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(100)# (period)
app.add(sort_method)

comp_info = gala.ComputeInfo(all_info, group)
DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')
DInfo.setPeriod(200)# (period)
app.add(DInfo)

xml = gala.XMLDump(all_info, 'p')
xml.setPeriod(100000)# (period)
xml.setOutputMass(True)
xml.setOutputImage(True)
xml.setOutputBond(True)
xml.setOutputAngle(True)
xml.setOutputInit(True)
xml.setOutputCris(True)
app.add(xml)

# groupA = gala.ParticleSet(all_info,'A')
# ect = gala.ExternalCenterTorque(all_info, groupA, 100.0)
# ect.setFieldDirection(1.0, 0.0, 0.0)
# ect.setPreNextShift(1, 2)
# app.add(ect)
#ready ro run 
app.run(10000)

Ebar=2.77
Ebind=6.0

# Period2: start the reaction
#DInfo.setPeriod(1)# (period)

rever = gala.DePolymerization(all_info, T, 16361)
rever.setParams('B-B', 1250.0, 1.0, 1.0, 1000.0, 180.0, Ebar+Ebind, 1.0, gala.DePolyFunc.harmonic)
# sets bondname, K_bond, r_0, b_0, K_angle, theta_0, epsilon0, Pr, and function.
rever.setPeriod(2000)
rever.setCountUnbonds(100000)
# sets how many steps to react.
app.add(rever)

reaction = gala.Polymerization(all_info, neighbor_list, 2**(1.0/6.0), 123)
reaction.setEnergyBar(Ebar)
reaction.setPr("B", "B", 1.0)
reaction.setNewBondTypeByPairs()
reaction.setNewAngleTypeByPairs()
reaction.generateAngle(True)
reaction.setMaxCris("B",1)
reaction.setInitInitReaction(True)
reaction.setAngleLowerLimitDegree(120)
reaction.setPeriod(2000)
app.add(reaction)

app.run(5000000)
# Period3: relaxation after the reaction
neighbor_list.printStats()
#
#
