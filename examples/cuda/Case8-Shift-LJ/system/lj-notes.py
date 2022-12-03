#!/usr/bin/python       #python 脚本声明
import cu_gala as gala 
from optparse import OptionParser      #导入option parser模块
 
global _options
parser = OptionParser()
parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
(_options, args) = parser.parse_args() #这四行是解析 --gpu=命令， 得到GPU的id

filename = 'relax.xml'    
build_method = gala.XMLReader(filename)  # 利用xml解析模块解读xml文件
perform_config = gala.PerformConfig(_options.gpu) #根据GPU ID来进行GPU绑定
all_info = gala.AllInfo(build_method, perform_config)  #建立体系信息
 
dt = 0.005
app = gala.Application(all_info, dt)  #建立应用

neighbor_list = gala.NeighborList(all_info, 4.0, 0.4)#(,rcut,rbuffer) #建立neigborlist，用来计算非键力
neighbor_list.addExclusionsFromBonds()   # 排除键连粒子的非键相互作用
neighbor_list.setFilterDiameters()       # 建立neiborlist的时候考虑粒子的diameter，也就是说截断半径是粒子的表面到表面的距离

slj0 = gala.SLJForce(all_info, neighbor_list, 4.0)# 建立shift LJ force计算，并设置参数
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

bondforce = gala.BondForceFENE(all_info) #建立键相作用，并设置参数
bondforce.setParams('A-C', 30.0, 1.5, 2.0, 1.0) #sets parameters by: bond type, K, r0, epsilon, sigma
bondforce.setParams('C-B', 30.0, 1.5, 2.0, 1.0)
bondforce.setParams('B-B', 30.0, 1.5, 2.0, 1.0)
bondforce.setConsiderDiameter(True)
app.add(bondforce)

group = gala.ParticleSet(all_info,'all') #设置粒子群组
comp_info = gala.ComputeInfo(all_info, group)#计算粒子群组的温度，压力等一些统计变量

bd=gala.BDNVT(all_info, group, 1.4, 123) # all_info, group, T, seed#建立布朗动力学积分方法
bd.setGamma('A', 3.0)
bd.setGamma('B', 1.0)
bd.setGamma('C', 1.0)
app.add(bd)

dinfo = gala.DumpInfo(all_info, comp_info, 'data.log') #输出体系的统计信息，温度，压力等等
dinfo.setPeriod(2000)
dinfo.dumpVirialEnergy(slj0)
dinfo.dumpVirialEnergy(slj1)
app.add(dinfo)

zm = gala.ZeroMomentum(all_info) #体系动量归零
zm.setPeriod(5000)
app.add(zm)

sort_method = gala.Sort(all_info) # 内存整理 ，用来加速计算
sort_method.setPeriod(200000)# (period)
app.add(sort_method)

mol2 = gala.MOL2Dump(all_info, 'particles') #输出mol2格式的构型文件
mol2.setPeriod(0)# (period)
app.add(mol2)

dcd = gala.DCDDump(all_info, 'particle',True)#输出DCD轨迹文件
dcd.setPeriod(200000)# (period)
dcd.unwrap(True)
app.add(dcd)

xml = gala.XMLDump(all_info, 'particle')#输出xml格式的构型文件
xml.setPeriod(200000)# (period)
xml.setOutputBond(True)
#xml.setOutputImage(True)
xml.setOutputDiameter(True)
xml.setOutputVelocity(True)
app.add(xml)

#ready ro run 
#app.run(10000)
#app.setDt(0.005)
app.run(300000001)#设置运行步数，并开始运行
neighbor_list.printStats()
#
#

