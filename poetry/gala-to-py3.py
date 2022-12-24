import os
import string
import math
import re
import platform

curr_path = os.path.abspath(__file__)
node = "/"
plat = platform.system().lower()
if plat == 'windows':
	print('Runing on windows')
	node = '\\'
p=curr_path.rfind(node)
if p!=-1:
	curr_path = curr_path[0:p]

texrp = [
['XmlReader', 'XMLReader'],
['Rnemd', 'RNEMD'],
['XmlDump', 'XMLDump'],
['Mol2Dump', 'MOL2Dump'],
['DcdDump', 'DCDDump'],
['BondForceFene', 'BondForceFENE'],
['DihedralForceOplsCosine', 'DihedralForceOPLSCosine'],
['SljForce', 'SLJForce'],
['LjConstrainForce', 'LJConstrainForce'],
['LjWallForce', 'LJWallForce'],
['LjForce', 'LJForce'],
['DpdForce', 'DPDForce'],
['DpdThermoLjForce', 'DPDThermoLJForce'],
['IntraMoleList', 'IntraMolList'],
['MdScfForce', 'MDSCFForce'],
['ItsMethod', 'ITSMethod'],
['DnaBuildXml', 'DNABuildXml'],
['DnaExForce', 'DNANoExForce'],
['DnaNoExForce', 'DNAExForce'],
['AniNpt', 'AniNPT'],
['BerendsenAniNvt', 'BerendsenAniNVT'],
['LzwForce', 'LZWForce'],
['NoseHooverAniNvt', 'NoseHooverAniNVT'],
['AndersenNvt', 'AndersenNVT'],
['BdNvt', 'BDNVT'],
['BerendsenNpt', 'BerendsenNPT'],
['BerendsenNvt', 'BerendsenNVT'],
['DpdGwvv', 'DPDGWVV'],
['LoweAndersenNvt', 'LoweAndersenNVT'],
['NoseHooverChainNvt', 'NoseHooverChainNVT'],
['NoseHooverNvt', 'NoseHooverNVT'],
['Npt', 'NPT'],
['Nve', 'NVE'],
['BdNvtRigid', 'BDNVTRigid'],
['BerendsenNptRigid', 'BerendsenNPTRigid'],
['NptRigid', 'NPTRigid'],
['NveRigid', 'NVERigid'],
['NvtRigid', 'NVTRigid'],
['galamost', 'gala'],
['addExclusionsFromBodys', 'addExclusionsFromBodies'],
]


for files in os.listdir(curr_path):
    if files.find(".gala") != -1:
        infile = files
        outfile = files[0:files.find(".gala")]+".py"
        print("convert ", infile," to " ,outfile)
        ingala = open(infile)
        outgala = open(outfile,"w")
        
        for line in ingala:
            if line.find("import sys") != -1:
                continue
            elif line.find("sys.path.append") != -1:
                continue
            elif line.find("OptionParser") != -1 or line.find("global _options") != -1 or line.find("parser.") != -1:
                continue
            elif line.find("import galamost") != -1:
                outgala.write("from poetry import cu_gala as gala\n")
                outgala.write("from poetry import _options\n")
                continue
                
            for j in texrp:
                strinfo = re.compile(j[0])
                line = strinfo.sub(j[1],line)
        #       print(line, j)
            outgala.write(line)
    