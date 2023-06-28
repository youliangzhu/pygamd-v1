import os
import string
import math
import re
import platform
# import logging

curr_path = os.getcwd()
'''
# 设置日志记录器
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

# 创建文件处理程序并添加到日志记录器
file_handler = logging.FileHandler('conversion.log')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# 输出当前路径
logger.info(f"Current Path: {os.getcwd()}")
'''

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
['BdNvt', 'LangevinNVT'],
['BerendsenNpt', 'BerendsenNPT'],
['BerendsenNvt', 'BerendsenNVT'],
['DpdGwvv', 'DPDGWVV'],
['LoweAndersenNvt', 'LoweAndersenNVT'],
['NoseHooverChainNvt', 'NoseHooverChainNVT'],
['NoseHooverNvt', 'NoseHooverNVT'],
['Npt', 'NPT'],
['Nve', 'NVE'],
['BdNvtRigid', 'LangevinNVTRigid'],
['BerendsenNptRigid', 'BerendsenNPTRigid'],
['NptRigid', 'NPTRigid'],
['NveRigid', 'NVERigid'],
['NvtRigid', 'NVTRigid'],
['galamost.', 'gala.'],
['addExclusionsFromBodys', 'addExclusionsFromBodies'],
['outPutXml', 'outPutXML'],
['outPutMol2', 'outPutMOL2'],
['Vsite\.VST', 'VST'],
['molgen\.Object.Shape', 'molgen.Shape'],
['molgen\.Protein.Model', 'molgen.Model'],
['molgen\.DNAchain.Strand', 'molgen.Strand'],
['PairForce\.Fun', 'PairFun'],
['Polymerization\.Func', 'PolyFunc'],
['DePolymerization\.Func', 'DePolyFunc'],
['PerformConfig\(int\(_options.gpu\)\)', 'PerformConfig(_options.gpu)'],
['PerformConfig\( int\(_options.gpu\)\)', 'PerformConfig(_options.gpu)'],
]

texrp = [[re.compile(pattern, re.IGNORECASE), replacement] for pattern, replacement in texrp]
'''
# 遍历当前目录中的文件
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
            # 替换文本
            for j in texrp:
                strinfo = re.compile(j[0])
                replaced_line = strinfo.sub(j[1], line)
                if replaced_line != line:
                    logger.info(f"Replaced '{j[0]}' with '{j[1]}': {replaced_line}")
                line = replaced_line
            outgala.write(line)
        ingala.close()
        outgala.close()
    elif files.find(".molg") != -1:
        infile = files
        outfile = files[0:files.find(".molg")]+"_new.molg"
        print("convert ", infile," to " ,outfile)
        ingala = open(infile)
        outgala = open(outfile,"w")
        
        for line in ingala:
            if line.find("import sys") != -1:
                continue
            elif line.find("sys.path.append") != -1:
                continue
            elif line.find("import molgen") != -1:
                outgala.write("from poetry import molgen\n")
                continue
                
            for j in texrp:
                strinfo = re.compile(j[0])
                line = strinfo.sub(j[1],line)
        #       print(line, j)
            outgala.write(line)
'''
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
    elif files.find(".molg") != -1:
        infile = files
        outfile = files[0:files.find(".molg")]+"_new.molg"
        print("convert ", infile," to " ,outfile)
        ingala = open(infile)
        outgala = open(outfile,"w")
        
        for line in ingala:
            if line.find("import sys") != -1:
                continue
            elif line.find("sys.path.append") != -1:
                continue
            elif line.find("import molgen") != -1:
                outgala.write("from poetry import molgen\n")
                continue
                
            for j in texrp:
                strinfo = re.compile(j[0])
                line = strinfo.sub(j[1],line)
        #       print(line, j)
            outgala.write(line)
