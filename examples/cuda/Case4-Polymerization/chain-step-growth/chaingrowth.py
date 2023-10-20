#!/usr/bin/python
from poetry import cu_gala as gala
from poetry import _options
import math


filename = 'gencon.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.002
T=1.0
app = gala.Application(all_info, dt)

all_info.addParticleType("D")
neighbor_list = gala.NeighborList(all_info, 2**(1.0/6.0), 0.2)#(,rcut,rbuffer) 
DPDThermoLJForce = gala.DPDThermoLJForce(all_info, neighbor_list, 2**(1.0/6.0), 12371)#(,,rcut,seed)
DPDThermoLJForce.setSigma(3.0)
DPDThermoLJForce.setParams('A', 'A', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'B', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('C', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('D', 'D', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'B', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'D', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'D', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('C', 'D', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
app.add(DPDThermoLJForce)


#Set bond potential(fene or harmonic) 
all_info.addBondType('sticky') 
bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('sticky', 1000.0 , 1.0)#(,K0, R0) # force constant, equlibrium position
app.add(bondforce)

all_info.addAngleType('sticky') 
angleforce = gala.AngleForceHarmonic(all_info)
angleforce.setParams('sticky', 10.0, 180.0)#(,K0, R0) # force constant, equlibrium position
app.add(angleforce)

# set temperature group
group=gala.ParticleSet(all_info, ['A', 'B', 'C', 'D'])

#computer information
comp_info = gala.ComputeInfo(all_info, group)

# integrate method
#nvt=gala.NoseHooverNVT(all_info, group,comp_info,T,tau) 
#app.add(nvt)
# NPT=gala.NPT(all_info, group,comp_info,comp_info,T,P,tau,tauP)
# app.add(NPT)
gwvv = gala.DPDGWVV(all_info, group)
app.add(gwvv)

#sort memory
sort_method = gala.Sort(all_info)
sort_method.setPeriod(300)# (period)
app.add(sort_method)

# ZeroMomentum
ZeroMomentum = gala.ZeroMomentum(all_info)
ZeroMomentum.setPeriod(100000)# (period)
app.add(ZeroMomentum)

DInfo = gala.DumpInfo(all_info, comp_info, 'data.log') #some basic infomation,temperature\momentum\pressure
DInfo.setPeriod(500)# (period)
app.add(DInfo)

# write information
xml = gala.XMLDump(all_info, 'p')
xml.setPeriod(100000)
xml.setOutputImage(True)
xml.setOutputBond(True)
xml.setOutputInit(True)
xml.setOutputCris(True)
xml.setOutputVelocity(True)
app.add(xml)


reaction1 = gala.Polymerization(all_info, neighbor_list, 2**(1.0/6.0) ,16361) 
reaction1.setPr('B', 'C', 0.5)
reaction1.setChangeTypeInReaction("B", "D")
reaction1.setChangeTypeInReaction("C", "B")
reaction1.setNewBondType("sticky")
reaction1.setNewAngleType("sticky")
reaction1.generateAngle(True)
reaction1.setPeriod(100)
app.add(reaction1)
xml.setPeriod(200)
app.run(100000)
neighbor_list.printStats()
