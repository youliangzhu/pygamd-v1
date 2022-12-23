#!/usr/bin/python
import os
from poetry import cu_gala as gala 
from poetry import _options

filename = 'topology.xml'# initial configuration file
build_method = gala.XMLReader(filename) #
perform_config = gala.PerformConfig(0)# assign GPU by index
all_info = gala.AllInfo(build_method, perform_config)# build system information

dt = 0.002  # timestep in ps
rcut = 1.5  # real-space cutoff in nm
rbuffer =0.1 # buffer-size (nm) for neighbour list
k_pme_mesh = 48
exclusionsSR='exclusionsSR.inc.py' # python script defining short-range exclusions
exclusionsEL='exclusionsEL.inc.py' # python script defining electrostatic exclusions
app = gala.Application(all_info, dt)# build up an application with system information and integration time-step
npoints_in_table = 1501 # Number of points in potential table

if npoints_in_table is None:
    raise ValueError('Value for number of points in the tabulated potential is not stated!!!')
# execute python scripts with exclusion definitions
exec(compile(source=open(exclusionsSR).read(), filename=exclusionsSR, mode='exec'))
exec(compile(source=open(exclusionsEL).read(), filename=exclusionsEL, mode='exec'))
# execute python script with tablulated potentials settings
exec(compile(source=open('tables.inc.py').read(), filename='tables.inc.py', mode='exec'))

group = gala.ParticleSet(all_info, "all")# Collection of all atoms in the system
group_charge = gala.ParticleSet(all_info, "charge") # Collection of charged atoms
#pppm = gala.PPPMForce(all_info, neighbor_listEL, group_charge)  # Turn on reciprocal part of PPPM
#pppm.setParams(k_pme_mesh, k_pme_mesh, k_pme_mesh, 5, rcut) # PPPM mesh 4x4x4, interpolate with 5th order polynom, real-space cutoff rcut
#app.add(pppm)
#kappa = pppm.getKappa() 
#print('kappa=', kappa)
#ewald = gala.EwaldForce(all_info, neighbor_listEL, group_charge, rcut) #Turn on real-space part of PPPM
#ewald.setParams(kappa)
#app.add(ewald)


comp_info = gala.ComputeInfo(all_info, group)  # calculating system informations, such as temperature, pressure, and momentum
Temperature = 303.0*0.00831  #kT  k=0.00831 kJ/mol

Bd = gala.BDNVT(all_info, group, Temperature, 123)# Set Langevin dynamics
Bd.setGamma(5.0)#(,gamma)
app.add(Bd)

ZeroMomentum = gala.ZeroMomentum(all_info) # removing the momentum of the center of mass
ZeroMomentum.setPeriod(100)# (period)
app.add(ZeroMomentum)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log')# output system informations, such as temperature, pressure, and momentum
DInfo.setPeriod(1000)# (period)
app.add(DInfo)

dcd = gala.DCDDump(all_info, 'trj',True) # write trajectory in dcd format
dcd.setPeriod(500)# (how often to write frames) 1ps
dcd.unwrap(True)
app.add(dcd)

#ready ro run 
app.run (5000000)#(the number of steps to run)
neighbor_listSR.printStats()# output the information about neighbor_list 
