#!/usr/bin/python
from poetry import cu_gala as gala 
from poetry import _options

filename = 'pn2d.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method, perform_config)

dt = 0.002
app = gala.Application(all_info, dt)

neighbor_list = gala.NeighborList(all_info, 1.0, 0.1)#(,rcut,rbuffer)
neighbor_list.addExclusionsFromBonds()
neighbor_list.addExclusionsFromAngles()

pf = gala.PairForce(all_info, neighbor_list)
pf.setParams('A', 'A',  25.0,  0.0, 0.0, 1.0, gala.PairForce.Func.harmonic) # A*epsilon, sigma, n,rcut
pf.setParams('B', 'B',  25.0,  0.0, 0.0, 1.0, gala.PairForce.Func.harmonic) # A*epsilon, sigma, n,rcut
pf.setParams('A', 'B',  25.0,  0.0, 0.0, 1.0, gala.PairForce.Func.harmonic) # A*epsilon, sigma, n,rcut
app.add(pf)

bf = gala.BondForceHarmonic(all_info)
bf.setParams('A-B', 1250.0, 1.0)#(,K0, R0)
#bf.setParams('B-B', 1250.0, 1.0)#(,K0, R0)
app.add(bf)

af = gala.AngleForceHarmonic(all_info)
af.setParams('B-A-B', 500.0, 120)#(,K0, R0)
#af.setParams('A-B-B', 500.0, 180)#(,K0, R0)
app.add(af)

group = gala.ParticleSet(all_info,'all')
comp_info = gala.ComputeInfo(all_info, group)


T = 1.0                                    #reduced unit
bd=gala.BDNVT(all_info, group, T, 123) # all_info, group, T, seed
app.add(bd)

sort_method = gala.Sort(all_info)
sort_method.setPeriod(100)# (period)
app.add(sort_method)

ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(100000)# (period)
app.add(ZeroMomentum)

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
# xml.setOutputInit(True)
# xml.setOutputCris(True)
app.add(xml)

#ready ro run 
app.run(50000000)
# Period3: relaxation after the reaction
neighbor_list.printStats()
#
#
