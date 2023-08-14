#!/usr/bin/python
from poetry import cu_gala as gala
from poetry import _options
from poetry import numerical
from poetry import molgen


# generate initial configuration
mol=molgen.Molecule(1)#number
mol.setParticleTypes("A")#type

gen=molgen.Generators(20, 20, 20)
gen.addMolecule(mol, 5000)#molecule,number
gen.setMinimumDistance(1.0)
gen.outPutXML("coor")

filename = 'coor.xml'# initial configuration file
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)# assign GPU by index
all_info = gala.AllInfo(build_method, perform_config)# build system information

dt = 0.001
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step
neighbor_list = gala.NeighborList(all_info, 3.0, 0.3)#(,rcut,rbuffer)

def lj(r, epsilon, sigma):
	v = 4.0 * epsilon * ( (sigma / r)**12 - (sigma / r)**6)
	return v

epsilon0 = 1.0
sigma0 = 1.0

pair = gala.PairForceTable(all_info, neighbor_list,  2000) # (,,the number of data points)
pair.setPotential('A', 'A' , numerical.pair(width=2000, func=lj, rmin=0.3, rmax=3.0, coeff=dict(epsilon=epsilon0, sigma=sigma0))) #(type,type,rcut,potential data file,starting column,ending column)
app.add(pair)

group = gala.ParticleSet(all_info, "all")# a collection of particles
comp_info = gala.ComputeInfo(all_info, group)  # calculating system informations, such as temperature, pressure, and momentum
T=1.0
P=0.1
NPT = gala.NPT(all_info, group, comp_info, comp_info, T, P, 1.0, 1.0)# temperature,pressure,tauT,tauP
app.add(NPT)

ZeroMomentum = gala.ZeroMomentum(all_info) # removing the momentum of the center of mass
ZeroMomentum.setPeriod(20000)# (period)
app.add(ZeroMomentum)

sort_method = gala.Sort(all_info) # memory sorting to improve data reading performance 
sort_method.setPeriod(1000)# (period)
app.add(sort_method)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')# output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(1000)# (period)
app.add(DInfo)

xml = gala.XMLDump(all_info, 'particle')# output the configuration files in xml format
xml.setPeriod(10000)# (period)
app.add(xml)

#ready ro run 
app.run (100000)#(the number of steps to run)
neighbor_list.printStats()# output the information about neighbor_list 
#
#
