from poetry import hip_gala as gala
import math

def exist(line ,name):
	if line.find(name)>=0:
		return True
	return False

nonbonded_params=[]
def parseNonbondParams(filename, atom_type):
	read_nonbond_params=False
	nonbond_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_nonbond_params and exist(nospace, '</pair_params>'):
			read_nonbond_params=False
			
		if read_nonbond_params:
			nonbond_params_file.append(la)

		if exist(nospace, '<pair_params>'):
			read_nonbond_params=True
	gala.close()

	for i in range(len(atom_type)):
		for j in range(i,len(atom_type)):
			for k in nonbond_params_file:
				if (k[0] == atom_type[i] and k[1] == atom_type[j]) or (k[0] == atom_type[j] and k[1] == atom_type[i]):
					nonbonded_params.append(k)			

bond_params = []
def parseBondParams(filename, bond_type):
	read_bond_params=False
	bond_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_bond_params and exist(nospace, '</bond_params>'):
			read_bond_params=False
			
		if read_bond_params:
			bond_params_file.append(la)

		if exist(nospace, '<bond_params>'):
			read_bond_params=True
	gala.close()

	for i in range(len(bond_type)):
		for k in bond_params_file:
			if (k[0] == bond_type[i]): 
				bond_params.append(k)	
	
angle_params = []
def parseAngleParams(filename, angle_type):
	read_angle_params=False
	angle_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_angle_params and exist(nospace, '</angle_params>'):
			read_angle_params=False
			
		if read_angle_params:
			angle_params_file.append(la)

		if exist(nospace, '<angle_params>'):
			read_angle_params=True
	gala.close()

	for i in range(len(angle_type)):
		for k in angle_params_file:
			if (k[0] == angle_type[i]):
				angle_params.append(k)	
	
dihedral_params = []
def parseDihedralParams(filename, dihedral_type):
	read_dihedral_params=False
	dihedral_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_dihedral_params and exist(nospace, '</dihedral_params>'):
			read_dihedral_params=False

		if read_dihedral_params:
			dihedral_params_file.append(la)

		if exist(nospace, '<dihedral_params>'):
			read_dihedral_params=True
	gala.close()

	for i in range(len(dihedral_type)):
		for k in dihedral_params_file:
			if (k[0] == dihedral_type[i]):
				dihedral_params.append(k)	
				
constraint_params = []
def parseConstraintParams(filename, constraint_type):
	read_constraint_params=False
	constraint_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_constraint_params and exist(nospace, '</constraint_params>'):
			read_constraint_params=False
			
		if read_constraint_params:
			constraint_params_file.append(la)

		if exist(nospace, '<constraint_params>'):
			read_constraint_params=True
	gala.close()

	for i in range(len(constraint_type)):
		for k in constraint_params_file:
			if (k[0] == constraint_type[i]):
				constraint_params.append(k)	

vsite_params = []
def parseVsiteParams(filename, vsite_type):
	read_vsite_params=False
	vsite_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_vsite_params and exist(nospace, '</vsite_params>'):
			read_vsite_params=False
			
		if read_vsite_params:
			vsite_params_file.append(la)

		if exist(nospace, '<vsite_params>'):
			read_vsite_params=True
	gala.close()

	for i in range(len(vsite_type)):
		for k in vsite_params_file:
			if (k[0] == vsite_type[i]):
				vsite_params.append(k)	

ah_params = []
def parseAHParams(filename, atom_type):
	read_ah_params=False
	ah_params_file=[]
	ah_params_file1=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_ah_params and exist(nospace, '</ah_params>'):
			read_ah_params=False
			
		if read_ah_params:
			ah_params_file.append(la)

		if exist(nospace, '<ah_params>'):
			read_ah_params=True
	gala.close()

	for i in range(len(atom_type)):
		for k in ah_params_file:
			if (k[0] == atom_type[i]):
				ah_params_file1.append(k)					

	for i in range(len(atom_type)):
		ahi = ah_params_file1[i]
		for j in range(i,len(atom_type)):
			ahj = ah_params_file1[j]
			ah_params.append([ahi[0], ahj[0], (float(ahi[1])+float(ahj[1]))/2.0,  (float(ahi[2])+float(ahj[2]))/2.0])	


wf_params = []
def parseWFParams(filename, atom_type):
	read_wf_params=False
	wf_params_file=[]
	gala=open(filename)
	for line in gala:		
		line1=""
		for cha in line:
			if cha!="#" and cha!=";":
				line1 += cha
			else:
				break
	
		lin=line1.strip('\n')
		la=lin.split()
		nospace=''.join(la)
		if len(la)==0 or la[0]==";":
			continue
			
		if read_wf_params and exist(nospace, '</wf_params>'):
			read_wf_params=False
			
		if read_wf_params:
			wf_params_file.append(la)

		if exist(nospace, '<wf_params>'):
			read_wf_params=True
	gala.close()

	for i in range(len(atom_type)):
		for j in range(i, len(atom_type)):    
			for k in wf_params_file:
				if (k[0] == atom_type[i] and k[1] == atom_type[j]) or (k[1] == atom_type[i] and k[0] == atom_type[j]) :
					wf_params.append(k)


def HarmonicForce(all_info, neighbor_list, alpha, filename):
	hf = gala.HarmonicForce(all_info, neighbor_list, 1.0)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(nonbonded_params)):
		hf.setParams(nonbonded_params[i][0], nonbonded_params[i][1], alpha, float(nonbonded_params[i][3])*1.2)
#		print pair[i][0], pair[i][1], value[i][0], value[i][1]
	return hf

def LJEwaldForce(all_info, neighbor_list, rcut, filename):
	lj = gala.LJEwaldForce(all_info, neighbor_list, rcut)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(nonbonded_params)):
		lj.setParams(nonbonded_params[i][0], nonbonded_params[i][1], float(nonbonded_params[i][2]), float(nonbonded_params[i][3]), float(nonbonded_params[i][4]))
	return lj
	
def LJCoulombShiftForce(all_info, neighbor_list, rcut, rshift, epsilonr, filename):
	lj = gala.LJCoulombShiftForce(all_info, neighbor_list)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	lj.setCoulomb(rcut, rshift, epsilonr)
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(nonbonded_params)):
		lj.setParams(nonbonded_params[i][0], nonbonded_params[i][1], float(nonbonded_params[i][2]), float(nonbonded_params[i][3]), float(nonbonded_params[i][4]), rcut, rshift)
	return lj	
	
	
def AHDHForce(all_info, neighbor_list, rcut, epsilon, debye_length, filename):
	ahdh = gala.AHDHForce(all_info, neighbor_list, rcut)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	ahdh.setDebyeLength(debye_length)
	parseAHParams(filename, atom_type)
	alpha = 1.0
	for i in range(0, len(ah_params)):
		ahdh.setParams(ah_params[i][0], ah_params[i][1], epsilon, float(ah_params[i][2]), alpha, float(ah_params[i][3]))
	return ahdh
    
def WFDHForce(all_info, neighbor_list, rcut, debye_length, filename):
	wfdh = gala.WFDHForce(all_info, neighbor_list, rcut)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	wfdh.setDebyeLength(debye_length)
	parseWFParams(filename, atom_type)
	for i in range(0, len(wf_params)):
		wfdh.setParams(wf_params[i][0], wf_params[i][1], float(wf_params[i][2]), float(wf_params[i][3]), float(wf_params[i][4]), float(wf_params[i][5]), float(wf_params[i][6]))
	return wfdh
	
def BondForceHarmonic(all_info, filename):
	bfh = gala.BondForceHarmonic(all_info)
	bond_type = all_info.getBondInfo().getBondTypes()	
	parseBondParams(filename, bond_type)
	for i in range(0, len(bond_params)):
#		print bond_params[i][0], float(bond_params[i][1]), float(bond_params[i][2])
		bfh.setParams(bond_params[i][0], float(bond_params[i][1]), float(bond_params[i][2]))
	return bfh

def AngleForceHarmonicCos(all_info, filename):
	afh = gala.AngleForceHarmonicCos(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_params)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0])	
		afh.setParams(angle_params[i][0], float(angle_params[i][1]), float(angle_params[i][2]))
	return afh
	
def AngleForceHarmonic(all_info, filename):
	afh = gala.AngleForceHarmonic(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_params)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0])	
		afh.setParams(angle_params[i][0], float(angle_params[i][1]), float(angle_params[i][2]))
	return afh	

def AngleForceUreyBradley(all_info, filename):
	afh = gala.AngleForceUreyBradley(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_params)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0])	
		afh.setParams(angle_params[i][0], float(angle_params[i][1]), float(angle_params[i][2]), float(angle_params[i][3]), float(angle_params[i][4]))
	return afh	

def DihedralForceAmberCosine(all_info, filename):
	dfh = gala.DihedralForceAmberCosine(all_info)
	dihedral_type = all_info.getDihedralInfo().getDihedralTypes()	
	parseDihedralParams(filename, dihedral_type)
	for i in range(0, len(dihedral_params)):
		dp = dihedral_params[i]
		if int(dp[9])==4:
#			print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), "gala.DihedralForceAmberCosine.Prop.improper"
			dfh.setParams(dp[0], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), gala.AmberProp.improper)
		elif int(dp[9])==9:
#			print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), "gala.DihedralForceAmberCosine.Prop.proper"
			dfh.setParams(dp[0], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), gala.AmberProp.proper) 
		else:
			raise RuntimeError('Error dihedral function type is not 4 or 9!')
	return dfh	
	
	
def DihedralForceRyckaertBellemans(all_info, filename):
	dfh = gala.DihedralForceRyckaertBellemans(all_info)
	dfh.setDividedFactorVDWELEC(0.5, 0.5)
	dihedral_type = all_info.getDihedralInfo().getDihedralTypes()	
	parseDihedralParams(filename, dihedral_type)
	for i in range(0, len(dihedral_params)):
		dp = dihedral_params[i]
#			print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), "gala.DihedralForceAmberCosine.Prop.improper"
		if int(dp[7])==3:
			dfh.setParams(dp[0], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]))
		else:
			raise RuntimeError('Error dihedral function type is not 3!')            
#		print dp
	return dfh		

def DihedralForceHarmonic(all_info, filename):
	dfh = gala.DihedralForceHarmonic(all_info)
	dihedral_type = all_info.getDihedralInfo().getDihedralTypes()	
	parseDihedralParams(filename, dihedral_type)
	for i in range(0, len(dihedral_params)):
		dp = dihedral_params[i]
		if int(dp[3])==2:
			angle=float(dp[1])
			if angle<0.0:
				angle += 360.0
			dfh.setParams(dp[0], float(dp[2]), angle, gala.HarmonicProp.improper)
		else:
			raise RuntimeError('Error dihedral function type is not 2!')
	return dfh		
	
def BondConstraint(all_info, filename):
	bc = gala.BondConstraint(all_info)
	constraint_type = all_info.getConstraintInfo().getConstraintTypes()	
	parseConstraintParams(filename, constraint_type)
	for i in range(0, len(constraint_params)):
		cp = constraint_params[i]
		bc.setParams(cp[0], float(cp[1]))
	return bc
	
def Vsite(all_info, filename):
	vs = gala.Vsite(all_info)
	vsite_type = all_info.getVsiteInfo().getVsiteTypes()	
	parseVsiteParams(filename, vsite_type)
	for i in range(0, len(vsite_params)):
		vp = vsite_params[i]
		if int(vp[4])==1:
			vs.setParams(vp[0], float(vp[1]), float(vp[2]), float(vp[3]), gala.VST.v3)
		elif int(vp[4])==4:
			vs.setParams(vp[0], float(vp[1]), float(vp[2]), float(vp[3]), gala.VST.v3out)
		else:
			raise RuntimeError('Error Vsite function type is not 1 and 4!')			
	return vs			


