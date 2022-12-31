from poetry import cu_gala as gala
import math
import string

param=[]
value=[]
pair=[]
define=[]
defaults=[]
atomproperties=[]

def exist(line ,name):
	if line.find(name)>=0:
		return True
	return False
	
def existType(name, atom_types):
	for i in atom_types:
		if i==name:
			return True
	return False

def get_sigma(name):
	for i in atomproperties:
		if i[0]==name:
			return float(i[len(i)-2])
def get_epsilon(name):
	for i in atomproperties:
		if i[0]==name:
			return float(i[len(i)-1])	

def parseNonbondParams(filename, atom_types):
	read_defaults=False
	read_nonbond_params=False
	read_pairtypes=False
	read_properties=False	
	itp=open(filename)
	for line in itp:
		if line.find("#define") != -1:
			lin=line.strip('\n')
			la=lin.split()
			if len(la)>=4:
				define.append((la[1],la[2],la[3]))
			
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
		if exist(nospace,'[defaults]'):
			read_defaults=True
			continue
		if read_defaults and len(la)>=2:
			defaults.append([la[0],la[1]])
			read_defaults=False	
		if read_nonbond_params and len(la)==5 and existType(la[0], atom_types) and existType(la[1], atom_types):
			pair.append((la[0],la[1]))
			param.append((la[3],la[4]))
		if read_nonbond_params and len(la)==4 and existType(la[0], atom_types) and existType(la[1], atom_types):
			findd=False
			for d in define:
				if d[0]==la[3]:
					pair.append((la[0],la[1]))
					param.append((d[1],d[2]))
					findd=True
					break
			if not findd:
				print (la[3])
				raise RuntimeError('Error interaction not found in define!')

		if read_pairtypes and len(la)==5 and existType(la[0], atom_types) and existType(la[1], atom_types):
			pair.append((la[0],la[1]))
			param.append((la[3],la[4]))
			# print "bbb", la
		if read_properties and len(la)>=6:
			atomproperties.append(la)
			
		if read_nonbond_params and len(la)<4:
			read_nonbond_params=False
		if read_pairtypes and len(la)<5:
			read_pairtypes=False
		if read_properties and len(la)<6:
			read_properties=False			
			
		if exist(nospace, '[nonbond_params]'):
			read_nonbond_params=True
		if exist(nospace, '[pairtypes]'):
			read_pairtypes=True
		if exist(nospace, '[atomtypes]'):
			read_properties = True			
	itp.close()
	if len(defaults)==0:
		defaults.append(("1","2"))
	
	if len(param)!=0:
		for i in range(0,len(param)):
			c6=string.atof(param[i][0])
			c12=string.atof(param[i][1])
			if c6!=0:
				sigma6 = c12/c6
				sigma = math.pow(sigma6,1.0/6.0) 
				epsilon4 = c6/sigma6
				epsilon = epsilon4/(4.0)
				if defaults[0][1]=='2':
					value.append((c12,c6))
				elif defaults[0][1]=='1':
					value.append((epsilon,sigma))
			else:
				value.append((0,0.47))
	# print value
	if len(param)==0:
		for i in range(len(atom_types)):
			for j in range(i,len(atom_types)):
				pair.append([atom_types[i],atom_types[j]])
				if defaults[0][1]=='2':
					sig=(get_sigma(atom_types[i])+get_sigma(atom_types[j]))/2
					eps=(get_epsilon(atom_types[i])*get_epsilon(atom_types[j]))**0.5
					value.append([eps,sig])	
	ffnb = open("nonboned-save.force_field","w")
	for i in range(len(pair)):
		ffnb.write(pair[i][0] + " "+ pair[i][1] +" "+ str(value[i][0]) +" "+ str(value[i][1]) + " 1.0\n")
	ffnb.close()

bond_params = []
def parseBondParams(filename, bond_types):
	read_bond_params = False
	bond_types_file = []
	itp=open(filename)
	for line in itp:
		if line.find("#define") != -1:
			lin=line.strip('\n')
			la=lin.split()
			if len(la)>=4:
				define.append((la[1],la[2],la[3]))
			
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
			
		if read_bond_params and len(la)<5:
			read_bond_params=False

		if read_bond_params:
			bond_types_file.append([la[0], la[1], la[2], la[3], la[4]])

		if exist(nospace, '[bondtypes]'):
			read_bond_params=True
	
	for bt in bond_types:
		la=bt.split("-")
		findd=False
		for btf in bond_types_file:
			if (btf[0]==la[0] and btf[1]==la[1]) or (btf[0]==la[1] and btf[1]==la[0]):
				findd=True
				bond_params.append([btf[3], btf[4], btf[2]])
				break
		if not findd:
			print (bt)
			raise RuntimeError('Error bond parameter not found in itp file!')
	itp.close()
	ffb = open("bond-save.force_field","w")
	for i in range(len(bond_types)):
		ffb.write(bond_types[i] +" "+ bond_params[i][1] +" "+ bond_params[i][0]+" "+ bond_params[i][2]+"\n")
	ffb.close()    
	
angle_params = []
def parseAngleParams(filename, angle_types):
	read_angle_params = False
	angle_types_file = []
	itp=open(filename)
	for line in itp:
		if line.find("#define") != -1:
			lin=line.strip('\n')
			la=lin.split()
			if len(la)>=4:
				define.append((la[1],la[2],la[3]))
			
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
			
		if read_angle_params and len(la)<6:
			read_angle_params=False

		if read_angle_params:
			if len(la)==6:
				angle_types_file.append([la[0], la[1], la[2], la[3], la[4], la[5]])
			elif len(la)==8:
				angle_types_file.append([la[0], la[1], la[2], la[3], la[4], la[5], la[6], la[7]])
			else:
				raise RuntimeError('Error angle data extract from itp file!')

		if exist(nospace, '[angletypes]'):
			read_angle_params=True
	
	for at in angle_types:
		la=at.split("-")
		findd=False
		for atf in angle_types_file:
			if (atf[0]==la[0] and atf[1]==la[1] and atf[2]==la[2]) or (atf[0]==la[2] and atf[1]==la[1] and atf[2]==la[0]):
				findd=True
				if len(atf)==6:
					angle_params.append([atf[4], atf[5], atf[3]])
				elif len(atf)==8:
					angle_params.append([atf[4], atf[5], atf[6], atf[7], atf[3]])
				break
		if not findd:
			print (at)
			raise RuntimeError('Error angle parameter not found in itp file!')
	itp.close()
    
	ffa = open("angle-save.force_field","w")
	for i in range(len(angle_types)):
		if len(angle_params[i]) == 3:
			ffa.write(angle_types[i] +" " + angle_params[i][1] +" "+ angle_params[i][0]+" "+ angle_params[i][2]+"\n")
		elif len(angle_params[i]) == 5:
			ffa.write(angle_types[i] +" " + angle_params[i][1] +" "+ angle_params[i][0]+" "+ angle_params[i][3]+" "+ angle_params[i][2]+" "+ angle_params[i][4]+"\n")            
	ffa.close()   	
	
dihedral_params = []
def parseDihedralParams(filename, dihedral_types):
	read_dihedral_params = False
	dihedral_types_file_read = []
	dihedral_types_file = []
	itp=open(filename)
	for line in itp:
		if line.find("#define") != -1:
			lin=line.strip('\n')
			la=lin.split()
			if len(la)>=4:
				define.append((la[1],la[2],la[3]))
			
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

		if read_dihedral_params and len(la)<8:
			read_dihedral_params=False

		if read_dihedral_params:
			dihedral_types_file_read.append(la)
			

		if exist(nospace, '[dihedraltypes]'):
			read_dihedral_params=True


	for k in range(0, len(dihedral_types_file_read)):
		if len(dihedral_types_file_read[k])==8 and dihedral_types_file_read[k][4] == "9":
			dtf = dihedral_types_file_read[k]
			dihedral_types_file_sort = dtf[:5]+['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0']
			if dtf[7]=='1':
				dihedral_types_file_sort[5]=dtf[6]
				dihedral_types_file_sort[9]=dtf[5]
			elif dtf[7]=='2':
				dihedral_types_file_sort[6]=dtf[6]
				dihedral_types_file_sort[10]=dtf[5]
			elif dtf[7]=='3':
				dihedral_types_file_sort[7]=dtf[6]
				dihedral_types_file_sort[11]=dtf[5]
			elif dtf[7]=='4':
				dihedral_types_file_sort[8]=dtf[6]
				dihedral_types_file_sort[12]=dtf[5]
			dihedral_types_file.append(dihedral_types_file_sort)
		else:
			dihedral_types_file.append(dihedral_types_file_read[k])


	for at in dihedral_types:
		al=at[:at.index("_F")]
		func=at[at.index("_F")+2]
		la=al.split("-")
		findd=False
#		print(al, la, func, at)
		for atf in dihedral_types_file:
			if ( (atf[0]==la[0] or atf[0]=='X') and (atf[1]==la[1] or atf[1]=='X')  and (atf[2]==la[2] or atf[2]=='X') and (atf[3]==la[3] or atf[3]=='X') and atf[4] == func) \
			or ( (atf[3]==la[0] or atf[3]=='X') and (atf[2]==la[1] or atf[2]=='X')  and (atf[1]==la[2] or atf[1]=='X') and (atf[0]==la[3] or atf[0]=='X') and atf[4] == func):
				findd=True
				dihedral_params.append(atf[4:]+[atf[4]])
				break
		if not findd:
			print (at)
			raise RuntimeError('Error dihedral parameter not found in itp file!')
	itp.close()	
	ffd = open("dihedral-save.force_field","w")
	for i in range(len(dihedral_types)):
		dp = ""
		for j in range(1, len(dihedral_params[i])):
			dp += dihedral_params[i][j]+" "
		ffd.write(dihedral_types[i] +" " + dp+ "\n")
	ffd.close()  

def LJEwaldForce(all_info, neighbor_list, rcut, filename):
	lj = gala.LJEwaldForce(all_info, neighbor_list, rcut)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(value)):
		lj.setParams(pair[i][0], pair[i][1], value[i][0], value[i][1], 1.0)
#		print pair[i][0], pair[i][1], value[i][0], value[i][1]
	return lj
    
def HarmonicForce(all_info, neighbor_list, alpha, filename):
	hf = gala.HarmonicForce(all_info, neighbor_list, 1.0)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(value)):
		hf.setParams(pair[i][0], pair[i][1], alpha, value[i][1]*1.2)
#		print pair[i][0], pair[i][1], value[i][0], value[i][1]
	return hf    
	
def LJCoulombShiftForce(all_info, neighbor_list, rcut, rshift, epsilonr, filename):
	lj = gala.LJCoulombShiftForce(all_info, neighbor_list)
	atom_type = all_info.getBasicInfo().getParticleTypes()	
	lj.setCoulomb(rcut, rshift, epsilonr)
	parseNonbondParams(filename, atom_type)
	for i in range(0, len(value)):
		lj.setParams(pair[i][0], pair[i][1], value[i][0], value[i][1], 1.0, rcut, rshift)
#		print pair[i][0], pair[i][1], value[i][0], value[i][1]
	return lj	
	
def BondForceHarmonic(all_info, filename):
	bfh = gala.BondForceHarmonic(all_info)
	bond_type = all_info.getBondInfo().getBondTypes()	
	parseBondParams(filename, bond_type)
	for i in range(0, len(bond_type)):
#		print bond_type[i], float(bond_params[i][1]), float(bond_params[i][0])
		bfh.setParams(bond_type[i], float(bond_params[i][1]), float(bond_params[i][0]))
	return bfh

def AngleForceHarmonicCos(all_info, filename):
	afh = gala.AngleForceHarmonicCos(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_type)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0])	
		afh.setParams(angle_type[i], float(angle_params[i][1]), float(angle_params[i][0]))
	return afh
	
def AngleForceHarmonic(all_info, filename):
	afh = gala.AngleForceHarmonic(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_type)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0])	
		afh.setParams(angle_type[i], float(angle_params[i][1]), float(angle_params[i][0]))
	return afh

def AngleForceUreyBradley(all_info, filename):
	afh = gala.AngleForceUreyBradley(all_info)
	angle_type = all_info.getAngleInfo().getAngleTypes()	
	parseAngleParams(filename, angle_type)
	for i in range(0, len(angle_type)):
#		print angle_type[i], float(angle_params[i][1]), float(angle_params[i][0]), float(angle_params[i][3]), float(angle_params[i][2])	
		afh.setParams(angle_type[i], float(angle_params[i][1]), float(angle_params[i][0]), float(angle_params[i][3]), float(angle_params[i][2]))
	return afh

def DihedralForceAmberCosine(all_info, filename):
	dfh = gala.DihedralForceAmberCosine(all_info)
	dihedral_type = all_info.getDihedralInfo().getDihedralTypes()	
	parseDihedralParams(filename, dihedral_type)
	for i in range(0, len(dihedral_type)):
		dp = dihedral_params[i]
		if int(dp[0])==4:
#			print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), "gala.AmberProp.improper"
			dfh.setParams(dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), gala.AmberProp.improper)
		elif int(dp[0])==9:
#			print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), "gala.DihedralForceAmberCosine.Prop.proper"
			dfh.setParams(dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]), float(dp[7]), float(dp[8]), gala.AmberProp.proper) 
		else:
			raise RuntimeError('Error dihedral function type is not 4 or 9!')
	return dfh	

def DihedralForceRyckaertBellemans(all_info, filename):
	dfh = gala.DihedralForceRyckaertBellemans(all_info)
	dfh.setDividedFactorVDWELEC(0.5, 0.5)
	dihedral_type = all_info.getDihedralInfo().getDihedralTypes()	
	parseDihedralParams(filename, dihedral_type)
	for i in range(0, len(dihedral_type)):
		dp = dihedral_params[i]
#		print dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6])), "gala.AmberProp.proper"
		dfh.setParams(dihedral_type[i], float(dp[1]), float(dp[2]), float(dp[3]), float(dp[4]), float(dp[5]), float(dp[6]))
	return dfh	
	


