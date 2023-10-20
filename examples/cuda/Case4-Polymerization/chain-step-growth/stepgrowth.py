#!/usr/bin/python
from poetry import cu_gala as gala
from poetry import _options
import math


filename = 'gencon.xml'
build_method = gala.XMLReader(filename)
perform_config = gala.PerformConfig(_options.gpu)
all_info = gala.AllInfo(build_method,perform_config)

dt = 0.002
#T=1.0;P=1.0;tau=0.5;tauP=1.0
app = gala.Application(all_info, dt)


all_info.addParticleType("C1")
all_info.addParticleType("B1")
neighbor_list = gala.NeighborList(all_info, 1.5, 0.2)#(,rcut,rbuffer) 
DPDThermoLJForce = gala.DPDThermoLJForce(all_info, neighbor_list, 1.12246 ,12371)#(,,rcut,seed)
DPDThermoLJForce.setSigma(3.0)
DPDThermoLJForce.setParams('A', 'A', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'B', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('C', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('C1', 'C1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B1', 'B1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)

DPDThermoLJForce.setParams('A', 'B', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'B1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('A', 'C1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)

DPDThermoLJForce.setParams('B', 'C', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'C1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('B', 'B1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)

DPDThermoLJForce.setParams('C', 'C1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
DPDThermoLJForce.setParams('C', 'B1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)

DPDThermoLJForce.setParams('C1', 'B1', 1.0 , 1.0 ,1.0)#(,,epsilon,sigma,alpha)
app.add(DPDThermoLJForce)


#Set bond potential(fene or harmonic) 
all_info.addBondType('sticky') 
bondforce = gala.BondForceHarmonic(all_info)
bondforce.setParams('sticky', 1000.0 , 1.0)#(,K0, R0) # force constant, equlibrium position
bondforce.setParams('B-B', 1000.0 , 1.0)#(,K0, R0) # force constant, equlibrium position
bondforce.setParams('C-C', 1000.0 , 1.0)#(,K0, R0) # force constant, equlibrium position
app.add(bondforce)
all_info.addAngleType('sticky') 
angleforce = gala.AngleForceHarmonic(all_info)
angleforce.setParams('sticky', 10.0, 180.0)#(,K0, R0) # force constant, equlibrium position
app.add(angleforce)

# set temperature group
group=gala.ParticleSet(all_info, ['A', 'B', 'C', 'C1', 'B1'])

#computer information
comp_info = gala.ComputeInfo(all_info, group)

# integrate method
#Gwvv = gala.DPDGWVV(all_info, group)
#app.add(Gwvv)
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

#write bin file
binary2 = gala.BinaryDump(all_info, 'particle')
binary2.setPeriod(10000)# (period)
binary2.setOutputForRestart()
app.add(binary2)

# write information
xml = gala.XMLDump(all_info, 'p')
xml.setPeriod(100000)
xml.setOutputImage(True)
xml.setOutputBond(True)
xml.setOutputInit(True)
xml.setOutputCris(True)
xml.setOutputVelocity(True)
app.add(xml)


# ljw=gala.LJWallForce(all_info, 0.8*2.0**(1.0/6.0))
# ljw.addWall(15.0, 0.0, 0.0, 1.0, 0.0, 0.0)
# ljw.addWall(0.0, 15.0, 0.0, 0.0, 1.0, 0.0)
# ljw.setParams("A", 1.0, 0.8, 1.0)
# ljw.setParams("B", 1.0, 0.8, 1.0)
# ljw.setParams("C", 1.0, 0.8, 1.0)
# ljw.setParams("B1", 1.0, 0.8, 1.0)
# ljw.setParams("C1", 1.0, 0.8, 1.0)
# app.add(ljw)


#Period3: start the SIP reaction, first depoly then poly to simulate the exchange reaction 
#reaction = gala.DePolymerization(all_info, T,16361)  
#reaction.setParams('sticky', 1111.0, 1.5, 0.967, 10, 1.0, gala.DePolyFunc.harmonic)
#(bondname, K, r_0, b_0, epsilon0, Pr,function:FENE,harmonic)
#reaction.setPeriod(100) # how many steps to react
#app.add(reaction)
app.run(1000)
reaction1 = gala.Polymerization(all_info, neighbor_list, 2**(1.0/6.0), 16361)
#(func_rule, K, r_0, b_0, epsilon0, function)
reaction1.setPr('C' , 'B', 0.1)
reaction1.setPr('C' ,'B1', 0.1)
reaction1.setPr('C1', 'B', 0.1)
reaction1.setPr('C1','B1', 0.1)
reaction1.setNewBondType("sticky")
reaction1.setNewAngleType("sticky")
reaction1.setMaxCris("C",  1)
reaction1.setMaxCris("B",  1)
reaction1.setMaxCris("C1", 1)
reaction1.setMaxCris("B1", 1)
reaction1.generateAngle(True)
reaction1.setChangeTypeInReaction("C", "C1")
reaction1.setChangeTypeInReaction("B", "B1")
reaction1.setPeriod(100)
app.add(reaction1)
xml.setPeriod(400)
app.run(200000)
neighbor_list.printStats()
