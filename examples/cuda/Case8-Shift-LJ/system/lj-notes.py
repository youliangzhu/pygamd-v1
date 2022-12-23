#!/usr/bin/python       #python �ű�����
from poetry import cu_gala as gala 
from poetry import _options

filename = 'relax.xml'    
build_method = gala.XMLReader(filename)  # ����xml����ģ����xml�ļ�
perform_config = gala.PerformConfig(_options.gpu) #����GPU ID������GPU��
all_info = gala.AllInfo(build_method, perform_config)  #������ϵ��Ϣ
 
dt = 0.005
app = gala.Application(all_info, dt)  #����Ӧ��

neighbor_list = gala.NeighborList(all_info, 4.0, 0.4)#(,rcut,rbuffer) #����neigborlist����������Ǽ���
neighbor_list.addExclusionsFromBonds()   # �ų��������ӵķǼ��໥����
neighbor_list.setFilterDiameters()       # ����neiborlist��ʱ�������ӵ�diameter��Ҳ����˵�ضϰ뾶�����ӵı��浽����ľ���

slj0 = gala.SLJForce(all_info, neighbor_list, 4.0)# ����shift LJ force���㣬�����ò���
slj0.setParams("A", "A", 0.0, 1.0, 0.1, 2.0+2**(1.0/6.0)) # type1, type2, epsilon, sigma, alpha, rcut
slj0.setParams("A", "C", 1.0, 1.0, 0.8, 4.0)
slj0.setParams("A", "B", 1.0, 1.0, 0.8, 4.0)
slj0.setParams("B", "B", 1.0, 1.0, 0.8, 3.0)
slj0.setParams("B", "C", 1.0, 1.0, 0.8, 3.0)
slj0.setParams("C", "C", 1.0, 1.0, 0.8, 3.0)
app.add(slj0)

slj1 = gala.SLJForce(all_info, neighbor_list, 4.0)# force rcut
slj1.setParams("A", "A", 0.5, 1.0, 0.1, 2.0+2**(1.0/6.0)) # type1, type2, epsilon, sigma, alpha, rcut
slj1.setParams("A", "C", 0.0, 1.0, 0.8, 4.0)
slj1.setParams("A", "B", 0.0, 1.0, 0.8, 4.0)
slj1.setParams("B", "B", 0.0, 1.0, 0.8, 3.0)
slj1.setParams("B", "C", 0.0, 1.0, 0.8, 3.0)
slj1.setParams("C", "C", 0.0, 1.0, 0.8, 3.0)
app.add(slj1)


#slj = gala.SLJForce(all_info, neighbor_list, 3.0)# force rcut
#slj.setParams("A", "A", 1.0, 1.0, 0.1, 1+2**(1.0/6.0)) # type1, type2, epsilon, sigma, alpha, rcut
#slj.setParams("A", "C", 1.0, 1.0, 0.0, 4.0)
#slj.setParams("A", "B", 0.5, 1.0, 1.0, 3.5)
#slj.setParams("C", "C", 1.0, 1.0, 0.8, 3.0)
#slj.setParams("C", "B", 1.0, 1.0, 0.8, 3.0)
#slj.setParams("B", "B", 1.0, 1.0, 1.0, 3.0)
#app.add(slj)

bondforce = gala.BondForceFENE(all_info) #�����������ã������ò���
bondforce.setParams('A-C', 30.0, 1.5, 2.0, 1.0) #sets parameters by: bond type, K, r0, epsilon, sigma
bondforce.setParams('C-B', 30.0, 1.5, 2.0, 1.0)
bondforce.setParams('B-B', 30.0, 1.5, 2.0, 1.0)
bondforce.setConsiderDiameter(True)
app.add(bondforce)

group = gala.ParticleSet(all_info,'all') #��������Ⱥ��
comp_info = gala.ComputeInfo(all_info, group)#��������Ⱥ����¶ȣ�ѹ����һЩͳ�Ʊ���

bd=gala.BDNVT(all_info, group, 1.4, 123) # all_info, group, T, seed#�������ʶ���ѧ���ַ���
bd.setGamma('A', 3.0)
bd.setGamma('B', 1.0)
bd.setGamma('C', 1.0)
app.add(bd)

dinfo = gala.DumpInfo(all_info, comp_info, 'data.log') #�����ϵ��ͳ����Ϣ���¶ȣ�ѹ���ȵ�
dinfo.setPeriod(2000)
dinfo.dumpVirialEnergy(slj0)
dinfo.dumpVirialEnergy(slj1)
app.add(dinfo)

zm = gala.ZeroMomentum(all_info) #��ϵ��������
zm.setPeriod(5000)
app.add(zm)

sort_method = gala.Sort(all_info) # �ڴ����� ���������ټ���
sort_method.setPeriod(200000)# (period)
app.add(sort_method)

mol2 = gala.MOL2Dump(all_info, 'particles') #���mol2��ʽ�Ĺ����ļ�
mol2.setPeriod(0)# (period)
app.add(mol2)

dcd = gala.DCDDump(all_info, 'particle',True)#���DCD�켣�ļ�
dcd.setPeriod(200000)# (period)
dcd.unwrap(True)
app.add(dcd)

xml = gala.XMLDump(all_info, 'particle')#���xml��ʽ�Ĺ����ļ�
xml.setPeriod(200000)# (period)
xml.setOutputBond(True)
#xml.setOutputImage(True)
xml.setOutputDiameter(True)
xml.setOutputVelocity(True)
app.add(xml)

#ready ro run 
#app.run(10000)
#app.setDt(0.005)
app.run(300000001)#�������в���������ʼ����
neighbor_list.printStats()
#
#

