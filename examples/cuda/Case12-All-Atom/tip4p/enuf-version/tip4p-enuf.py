#!/usr/bin/python
import cu_gala as gala 
from optparse import OptionParser

global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args()

filename = 'tip4p.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.001
app = gala.Application(all_info, dt)
neighbor_list = gala.NeighborList(all_info, 0.9, 0.09)#(,rcut,rbuffer)
neighbor_list.exclusion(["constraint", "vsite"])

lj = gala.LJEwaldForce(all_info, neighbor_list, 0.9)#(,,rcut)
lj.setParams('OW', 'OW',   0.648520, 0.315365, 1.0)
lj.setParams('HW', 'HW',   0.0, 0.47, 1.0)
lj.setParams('MW',  'MW',  0.0, 0.47, 1.0)

lj.setParams('OW', 'HW',   0.0, 0.47, 1.0)
lj.setParams('OW',  'MW',  0.0, 0.47, 1.0)

lj.setParams('HW',  'MW',  0.0, 0.47, 1.0)
lj.setEnergy_shift()
lj.setDispVirialCorr(True)#dispersion virial correction
app.add(lj) 

bondconstraint = gala.BondConstraint(all_info)# bond constraints using LINCS algorithm
bondconstraint.setParams('oh', 0.09572)#(type, r0)
bondconstraint.setParams('hh', 0.15139)#(type, r0)
bondconstraint.setExpansionOrder(4)
bondconstraint.setNumIters(2)
app.add(bondconstraint)

groupC = gala.ParticleSet(all_info, "charge")
kappa=3.05796313286
enuf = gala.ENUFForce(all_info, neighbor_list, groupC)
enuf.setParams(kappa, 2.0, 2, 32, 32, 32) #  alpha, sigma, precision, grid number in x, y, and z directions
app.add(enuf)

vs = gala.Vsite(all_info)#virtual interaction sites
vs.setParams('v', 0.128012065, 0.128012065, 0.0, gala.Vsite.VST.v3 )
app.add(vs)

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
DInfo.dumpVirialEnergy(lj)
DInfo.dumpVirialEnergy(enuf)
app.add(DInfo)

xml = gala.XMLDump(all_info, 'p')
xml.setPeriod(10000)# (period)
xml.setOutput(["constraint", "vsite", "mass", "velocity", "charge","image"])
app.add(xml)

# Period1
#ready ro run 
app.run(1000000)
neighbor_list.printStats()
#
#
