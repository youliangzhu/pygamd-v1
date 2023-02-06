import string
import math

from optparse import OptionParser
 
global _options
parser = OptionParser()
parser.add_option('--top', dest='top',help='top file')
parser.add_option('--gro', dest='gro',help='gro file')
parser.add_option('--er', dest='er',help='dielectric constant')
(_options, args) = parser.parse_args()


topfile=_options.top
grofile=_options.gro
er=_options.er

if not topfile:
    raise RuntimeError("Error, please set top file by '--top=' !")  
if not grofile:
    raise RuntimeError("Error, please set top file by '--gro=' !")      
if not er:
    er=1.0
else:
    er=float(er)    
print("dielectric constant = ", er)
#---------------
# read gro file
#---------------
ItpFiles=[]
Molecules=[]

def exist(line ,name):
    if line.find(name)>=0:
        return True
    return False
    
def existSeparated(array, name):
    for i in array:
        if i==name:
            return True
    return False

pos=[]
vel=[]
box=[]

def ParseGro(filename):
    print ("read" + filename)
    count=0 
    gro=open(filename)
    l0=gro.readline()
    l1=gro.readline()
    NP = int(l1)
    for l in range(0, NP):
        line = gro.readline()
        s1=line[0:15]
        s2=line[15:44]
        s3=line[44:52]
        s4=line[52:60]
        s5=line[60:]    
        line=s1+" "+s2+" "+s3+" "+s4+" "+s5
        lin=line.strip('\n')
#       print lin       
        la=lin.split()
        if len(la)==0 or la[0]==";":
            continue
        if len(la)==6:
            pos.append([float(la[3]), float(la[4]), float(la[5])])
            count += 1
        if len(la)==9:
            pos.append([float(la[3]), float(la[4]), float(la[5])])
            vel.append([float(la[6]), float(la[7]), float(la[8])])
            count += 1
    l2=gro.readline()
    lin=l2.strip('\n')
    la=lin.split()
    lx=float(la[0])
    ly=float(la[1])
    lz=float(la[2])
    box.append([lx, ly, lz])
    gro.close()
    if NP != count:
        raise RuntimeError('Error read line number is wrong!');
    
ParseGro(grofile)   
lx=box[0][0]
ly=box[0][1]
lz=box[0][2]


#---------------
# read top file
#---------------
def ParseTop(filename):
    print ("read" + filename)
    top=open(filename)
    for line in top:
        lin=line.strip('\n')
        la=lin.split()
        if len(la)==0 or la[0]==";":
            continue
        if len(la)>=2 and la[0]=="#include":
            n=len(la[1])
            ItpFiles.append(la[1][1:n-1])
    top.seek(0,0)
    read=False
    for line in top:
        lin=line.strip('\n')
        la=lin.split()
        nospace=''.join(la)
        if len(la)==0 or la[0]==";":
            continue        
        if read and len(la)==2:
            Molecules.append((la[0], int(la[1])))
        if exist(nospace, "[molecules]"):
            read=True
    top.close()
ParseTop(topfile)

#---------------------------------------------------------
#  read bondtypes, constrainttypes, angletypes, dihedraltypes
#---------------------------------------------------------

#---------------
# read atomproperty
#---------------
atomproperties=[]
bondtypes=[]
constrainttypes=[]
angletypes=[]
dihedraltypes=[]

read_atomproperties=False
read_bondtypes=False
read_constraintstypes=False
read_angletypes=False
read_dihedraltypes=False

for filename in ItpFiles:
    print("read itpfile ",filename)
    itp=open(filename,'r')
    for line in itp:
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
        if read_atomproperties and len(la)<6:
            read_atomproperties = False             
        if read_bondtypes and len(la)<5:
            read_bondtypes=False
        if read_constraintstypes and len(la)<4:
            read_constraintstypes =False
        if read_angletypes and len(la)<6:
            read_angletypes=False
        if read_dihedraltypes and len(la)<8:
            read_dihedraltypes=False            

            
        if read_atomproperties and len(la)>=6:
            atomproperties.append(la)            
        if read_bondtypes and len(la)==5:
            bondtypes.append(la)
            # print la
        if read_constraintstypes and len(la)==4:
            constrainttypes.append(la)
            # print la
        if read_angletypes and (len(la)==6 or len(la)==8):
            angletypes.append(la)
            # print la
        if read_dihedraltypes and (len(la)==8 or len(la)==11):
            dihedraltypes.append(la)
            # print la
        

        if exist(nospace, "[atomtypes]"):
            read_atomproperties = True
            read_bondtypes=False
            read_constraintstypes=False
            read_angletypes=False
            read_dihedraltypes=False            
        if  exist(nospace, "[bondtypes]"):
            read_atomproperties = False        
            read_bondtypes=True
            read_constraintstypes=False
            read_angletypes=False
            read_dihedraltypes=False
        if  exist(nospace, "[constrainttypes]"):
            read_atomproperties = False         
            read_bondtypes=False
            read_constraintstypes=True
            read_angletypes=False
            read_dihedraltypes=False
        if  exist(nospace, "[angletypes]"):
            read_atomproperties = False         
            read_bondtypes=False
            read_constraintstypes=False
            read_angletypes=True
            read_dihedraltypes=False
        if  exist(nospace, "[dihedraltypes]"):
            read_atomproperties = False         
            read_bondtypes=False
            read_constraintstypes=False
            read_angletypes=False
            read_dihedraltypes=True
    itp.close()

#------------------------------------------------------------------
#  sort dihedraltypes
#------------------------------------------------------------------

def find5(array,name):
    for i in array:
        if i[0]==name[0] and i[1]==name[1] and i[2]==name[2] and i[3]==name[3] and i[4]==name[4]:
            return True
    return False

dihedral_types=[]
for i in dihedraltypes:
    if not find5(dihedral_types,i[:5]):
        dihedral_types.append(i[:5]+['0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0'])

for i in range(len(dihedral_types)):
    for k in dihedraltypes:
        if dihedral_types[i][:5]==k[:5]:
            if k[4] == '9':
                if k[7]=='1':
                    dihedral_types[i][5]=k[6]
                    dihedral_types[i][9]=k[5]
                elif k[7]=='2':
                    dihedral_types[i][6]=k[6]
                    dihedral_types[i][10]=k[5]
                elif k[7]=='3':
                    dihedral_types[i][7]=k[6]
                    dihedral_types[i][11]=k[5]
                elif k[7]=='4':
                    dihedral_types[i][8]=k[6]
                    dihedral_types[i][12]=k[5]
            elif k[4] == '3':
                dihedral_types[i][5:]=k[5:]            
#print (dihedraltypes)
#print (dihedral_types)

def bond_parameter(atom1,atom2,type1):
    for i in bondtypes:
        if (i[0]==atom1 and i[1]==atom2 and i[2]==type1) or (i[1]==atom1 and i[0]==atom2 and i[2]==type1) :
            return i[2:5]
    return False

def angle_parameter(atom1,atom2,atom3,type1):
    for i in angletypes:
        if (i[0]==atom1 and i[1]==atom2 and i[2]==atom3 and i[3]==type1) or (i[2]==atom1 and i[1]==atom2 and i[0]==atom3 and i[3]==type1) :
            if len(i)==6:
                return i[3:6]
            elif len(i)==8:
                return i[3:8]
    return False

def dihedral_parameter(atom1,atom2,atom3,atom4,type1):
    for i in dihedral_types:
        if ( (i[0]==atom1 or i[0]=='X') and (i[1]==atom2 or i[1]=='X')  and (i[2]==atom3 or i[2]=='X') and (i[3]==atom4 or i[3]=='X') and i[4]==type1) \
        or ( (i[3]==atom1 or i[3]=='X') and (i[2]==atom2 or i[2]=='X')  and (i[1]==atom3 or i[1]=='X') and (i[0]==atom4 or i[0]=='X') and i[4]==type1):
            return i[4:13]
    return False


#------------------------
#parse molecules
#------------------------
nametotype = {}
if (len(atomproperties[0])==8):
    for i in atomproperties:
        nametotype[i[0]] = i[1]
    print("associate name to type")
#    print("associate name to type", nametotype)    
atom_all=[]
bond_all=[]
constraint_all=[]
pairs_all=[]
virtual_all=[]
angle_all=[]
dihedral_all=[]
exclusion_all=[]
alia = []
for i in Molecules:
    print (str(i))
    mol_name=i[0]
    mol_num=i[1]
    read_molecule=False
    read_atom=False
    read_bond=False
    read_constraints=False
    read_pair=False
    read_virtual=False
    read_angle=False
    read_dihedral=False
    read_exclusion=False
    atom=[]
    bond=[]
    constraint=[]
    pairs=[]
    virtual=[]  
    angle=[]
    dihedral=[]
    exclusion=[]
    for filename in ItpFiles:
        itp=open(filename)
        for line in itp:
            if line.find("#define") != -1 :
                lin=line.strip('\n')
                la=lin.split()
                if len(la)>=4:
                    alia.append((la[1],la[2],la[3]))
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
            if la[0]=='[':
                read_atom=False
                read_bond=False
                read_constraints=False
                read_pair=False
                read_virtual=False
                read_angle=False
                read_dihedral=False
                read_exclusion=False
            if read_atom and len(la)>=7:
                atom.append(la)
                # print la

            if read_bond and len(la)==4:
                findd=False
                for a in alia:
                    if a[0]==la[3]:
                        bond.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]]+la[:3]+[a[1],a[2]])
                        # print la[0], la[1], la[2], a[1], a[2]
                        findd=True
                        break
                if not findd:
                    # print la[3]
                    raise RuntimeError('Error bond not found in alia!') 
            elif read_bond and len(la)==3:
                if bool(bond_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],la[2])):
                    para=bond_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],la[2])
                    bond.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]]+la[:2]+para)
                elif bool(bond_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]],la[2])):
                    para=bond_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]],la[2])
                    bond.append([nametotype[atom[int(la[0])-1][1]]+"-"+nametotype[atom[int(la[1])-1][1]]]+la[:2]+para)
                else:
                    print(la)
                    print(atom[int(la[0])-1][1],atom[int(la[1])-1][1], la[2])
                    print(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], la[2]) 
                    raise RuntimeError('Error bond not found in forcefield!')   
            elif read_bond and len(la)>=5:
                bond.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]]+la)
            
            
            if read_constraints and len(la)>=4:
                constraint.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]]+la)
                # print(la, "constraint")

            if read_pair and len(la)>=3:
                pairs.append(la)
                # print(la)

            if read_virtual and len(la)>=7:
                if len(la)==7:
                    la.append('0.0')
                virtual.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]+"-"+atom[int(la[3])-1][1]]+la)
                # print(la)

            if read_angle and len(la)==5:
                findd=False
                for a in alia:
                    if a[0]==la[4]:
                        angle.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]] + la[:4] + [a[1], a[2]])
                        # print la[0], la[1], la[2], la[3], a[1], a[2]
                        findd=True
                        break
                if not findd:
                    # print la[4]
                    raise RuntimeError('Error bond not found in alia!')     
            elif read_angle and len(la)==4:     
                if bool(angle_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1] ,la[3])):
                    para=angle_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1] ,la[3])
                    angle.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]]+la[:3]+para)
                elif bool(angle_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], la[3])):
                    para=angle_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], la[3])
                    angle.append([nametotype[atom[int(la[0])-1][1]]+"-"+nametotype[atom[int(la[1])-1][1]]+"-"+nametotype[atom[int(la[2])-1][1]]]+la[:3]+para)                    
                else:
                    print(la)
                    print(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1], la[3])
                    print(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], la[3])             
                    raise RuntimeError('Error angle not found in forcefield!')
            elif read_angle and len(la)>=6: 
                angle.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]]+la)            
        

            if read_dihedral and len(la)>=7:
                dihedral.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]+"-"+atom[int(la[3])-1][1]+"_F"+la[4]]+la)
                # print la
            elif read_dihedral and len(la)==5:
                if bool(dihedral_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1] ,atom[int(la[3])-1][1],la[4])):
                    para=dihedral_parameter(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1] ,atom[int(la[3])-1][1],la[4])
                    dihedral.append([atom[int(la[0])-1][1]+"-"+atom[int(la[1])-1][1]+"-"+atom[int(la[2])-1][1]+"-"+atom[int(la[3])-1][1]+"_F"+la[4]]+la[:4]+para)
                elif bool(dihedral_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], nametotype[atom[int(la[3])-1][1]], la[4])):
                    para=dihedral_parameter(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], nametotype[atom[int(la[3])-1][1]], la[4])
                    dihedral.append([nametotype[atom[int(la[0])-1][1]]+"-"+nametotype[atom[int(la[1])-1][1]]+"-"+nametotype[atom[int(la[2])-1][1]]+"-"+nametotype[atom[int(la[3])-1][1]]+"_F"+la[4]]+la[:4]+para)                    
                else:
                    print(la)
                    print(atom[int(la[0])-1][1],atom[int(la[1])-1][1],atom[int(la[2])-1][1] ,atom[int(la[3])-1][1],la[4])
                    print(nametotype[atom[int(la[0])-1][1]], nametotype[atom[int(la[1])-1][1]], nametotype[atom[int(la[2])-1][1]], nametotype[atom[int(la[3])-1][1]], la[4])
                    raise RuntimeError('Error dihedral not found in forcefield!')
                    
            if read_exclusion:
                for pj in range(1, len(la)):
                    if int(la[0])>0 and int(la[pj]) >0:
                        exclusion.append([int(la[0])-1, int(la[pj])-1])
                    else:
                        raise RuntimeError('Error exclusion index!')
                    
                
            if read_atom and len(la)<7:
                read_atom = False
            if read_bond and len(la)<3:
                read_bond = False   
            if read_pair and len(la)<3:
                read_pair = False
            if read_constraints and len(la)<4:
                read_constraints = False    
            if read_virtual and len(la)<7 :
                read_virtual=False              
            if read_angle and len(la)<4:
                read_angle = False
            if read_dihedral and len(la)<5:
                read_dihedral = False
            if read_exclusion and len(la)<2:
                read_exclusion = False
            if read_molecule and exist(nospace, "[moleculetype]"):
                read_molecule = False
            
            if len(la)==2 and existSeparated(la, mol_name):
                read_molecule=True
            if read_molecule and exist(nospace, "[atoms]"):
                read_atom=True
            if read_molecule and exist(nospace, "[bonds]"):
                read_bond=True
                read_atom=False
            if read_molecule and exist(nospace, "[pairs]"):
                read_bond=False
                read_pair=True
            if read_molecule and exist(nospace, "[constraints]"):
                read_constraints=True
                read_atom=False
                read_bond=False 
            if read_molecule and exist(nospace,'[virtual_sites3]'):
                read_virtual=True
                read_constraints=False
                read_atom=False
                read_bond=False     
            if read_molecule and exist(nospace, "[angles]"):
                read_angle=True
                read_bond=False
                read_virtual=False
                read_constraints=False
                read_atom=False
            if read_molecule and exist(nospace, "[dihedrals]"):
                read_dihedral=True
                read_angle=False
            if read_molecule and exist(nospace, "[exclusions]"):
                read_exclusion=True
                read_bond=False
                read_angle=False
                read_virtual=False
                read_constraints=False
                read_dihedral=False

        itp.close()
    atom_all.append(atom)
    bond_all.append(bond)
    constraint_all.append(constraint)
    virtual_all.append(virtual) 
    angle_all.append(angle)
    dihedral_all.append(dihedral)
    exclusion_all.append(exclusion)
# print (alia)

#---------------
# sort bonds, angles and dihedrals
#---------------    
atom_all_new=[]
bond_all_new=[]
constraint_all_new=[]
virtual_all_new=[]
angle_all_new=[]
dihedral_all_new=[]


atom_type=[]
bond_type=[]
constraint_type=[]
virtual_type=[]
angle_type=[]
dihedral_type=[]


def existBondedType(array, name):
    for i in array:
        if i[0]==name:
            return True
    return False
    
def existStrictBondedType(array, name): 
    for i in array:
        if i[0].find(name[0])>=0 and i[1:]==name[1:]:
            return i[0]

    
def existAtomType(array, name):
    for i in array:
        if i==name:
            return True
    return False    


bond_map={}
constraint_map={}
virtual_map={}
angle_map={}
dihedral_map={} 
for i in range(0, len(Molecules)):
    mol_name=Molecules[i][0]
    mol_num=Molecules[i][1]
    atom =atom_all[i]
    bond = bond_all[i]
    constraint = constraint_all[i]  
    virtual = virtual_all[i]
    angle = angle_all[i]
    dihedral = dihedral_all[i]
    atom_new=[]
    bond_new=[]
    constraint_new=[]
    virtual_new=[]
    angle_new=[]  
    dihedral_new=[]
    for j in atom:
        atom_new.append(j[1])
        if not existAtomType(atom_type, j[1]):
            atom_type.append(j[1])
    for k in bond:
        if not existBondedType(bond_type, k[0]):
            bond_type.append([k[0], k[3], k[4], k[5]])
            bond_new.append([k[0], int(k[1])-1, int(k[2])-1])
            bond_map.update({k[0]:0})
        else:
            name = existStrictBondedType(bond_type, [k[0], k[3], k[4], k[5]])
            if name is None:
                bond_map[k[0]] += 1
                bond_type.append([k[0]+"_"+str(bond_map[k[0]]), k[3], k[4], k[5]])
                bond_new.append([k[0]+"_"+str(bond_map[k[0]]), int(k[1])-1, int(k[2])-1])
            else:
                bond_new.append([name, int(k[1])-1, int(k[2])-1])

    for o in constraint:
        if not existBondedType(constraint_type, o[0]):
            constraint_type.append([o[0], o[3], o[4]])
            constraint_new.append([o[0], int(o[1])-1, int(o[2])-1])
            constraint_map.update({o[0]:0})
        else:
            name = existStrictBondedType(constraint_type, [o[0], o[3], o[4]])
            if name is None:
                constraint_map[o[0]] += 1
                constraint_type.append([o[0]+"_"+str(constraint_map[o[0]]), o[3], o[4]])
                constraint_new.append([o[0]+"_"+str(constraint_map[o[0]]), int(o[1])-1, int(o[2])-1])
            else:
                constraint_new.append([name, int(o[1])-1, int(o[2])-1])     

    for p in virtual:
        if not existBondedType(virtual_type, p[0]):
            virtual_type.append([p[0], p[5], p[6], p[7], p[8]])
            virtual_new.append([p[0], int(p[1])-1, int(p[2])-1, int(p[3])-1,int(p[4])-1])
            virtual_map.update({p[0]:0})
        else:
            name = existStrictBondedType(virtual_type, [p[0], p[5], p[6], p[7], p[8]])
            if name is None:
                virtual_map[p[0]] += 1
                virtual_type.append([p[0]+"_"+str(virtual_map[p[0]]), p[5], p[6], p[7], p[8]])
                virtual_new.append([p[0]+"_"+str(virtual_map[p[0]]), int(p[1])-1, int(p[2])-1, int(p[3])-1,int(p[4])-1])
            else:
                virtual_new.append([name, int(p[1])-1, int(p[2])-1, int(p[3])-1,int(p[4])-1])       

    # for l in angle:
        # if not existBondedType(angle_type, l[0]):
            # angle_type.append([l[0], l[4], l[5], l[6]])
        # angle_new.append([l[0], int(l[1])-1, int(l[2])-1, int(l[3])-1])
        
    for l in angle:
        if not existBondedType(angle_type, l[0]):
            if len(l)==7:
                angle_type.append([l[0], l[4], l[5], l[6]])
            elif len(l)==9:
                angle_type.append([l[0], l[4], l[5], l[6], l[7], l[8]])         
            angle_new.append([l[0], int(l[1])-1, int(l[2])-1, int(l[3])-1])
            angle_map.update({l[0]:0})
        else:
            if len(l)==7:
                name = existStrictBondedType(angle_type, [l[0], l[4], l[5], l[6]])
                if name is None:
                    angle_map[l[0]] += 1            
                    angle_type.append([l[0]+"_"+str(angle_map[l[0]]), l[4], l[5], l[6]])            
                    angle_new.append([l[0]+"_"+str(angle_map[l[0]]), int(l[1])-1, int(l[2])-1, int(l[3])-1])
                else:
                    angle_new.append([name, int(l[1])-1, int(l[2])-1, int(l[3])-1])             
            elif len(l)==9:
                name = existStrictBondedType(angle_type, [l[0], l[4], l[5], l[6], l[7], l[8]])
                if name is None:
                    angle_map[l[0]] += 1
                    angle_type.append([l[0]+"_"+str(angle_map[l[0]]), l[4], l[5], l[6], l[7], l[8]])                
                    angle_new.append([l[0]+"_"+str(angle_map[l[0]]), int(l[1])-1, int(l[2])-1, int(l[3])-1])
                else:
                    angle_new.append([name, int(l[1])-1, int(l[2])-1, int(l[3])-1])             
    

    # for m in dihedral:
        # if not existBondedType(dihedral_type, m[0]):
            # dihedral_type.append([m[0]]+m[5:])
        # dihedral_new.append([m[0], int(m[1])-1, int(m[2])-1, int(m[3])-1, int(m[4])-1])

    for m in dihedral:
        if not existBondedType(dihedral_type, m[0]):
            dihedral_type.append([m[0]]+m[5:])
            dihedral_new.append([m[0], int(m[1])-1, int(m[2])-1, int(m[3])-1, int(m[4])-1])
            dihedral_map.update({m[0]:0})
        else:
            name = existStrictBondedType(dihedral_type, [m[0]]+m[5:])
            if name is None:
                dihedral_map[m[0]] += 1
                dihedral_type.append([m[0]+"_"+str(dihedral_map[m[0]])]+m[5:])
                dihedral_new.append([m[0]+"_"+str(dihedral_map[m[0]]), int(m[1])-1, int(m[2])-1, int(m[3])-1, int(m[4])-1])
            else:
                dihedral_new.append([name, int(m[1])-1, int(m[2])-1, int(m[3])-1, int(m[4])-1])         

    atom_all_new.append(atom_new)
    bond_all_new.append(bond_new)
    constraint_all_new.append(constraint_new)
    virtual_all_new.append(virtual_new) 
    angle_all_new.append(angle_new)
    dihedral_all_new.append(dihedral_new)

#---------------
# sort mass
#---------------
def findMass(name):
    for j in atomproperties:
        if j[0]==name:
            if len(j)==7:
                return j[2]
            elif len(j)==6:
                return j[1]
    print (str(name))
    raise RuntimeError('Error particle mass can not be found!')

mass_all_new=[]
for i in range(0, len(Molecules)):
    atom=atom_all[i]
    if len(atom[0])==7:
        atom =atom_all_new[i]
        mass=[]
        for j in atom:
            mass.append(findMass(j))
        mass_all_new.append(mass)
    elif len(atom[0])==8:
        mass=[]
        for j in range(0, len(atom)):
            mass.append(atom[j][7])
        mass_all_new.append(mass)
#---------------
# sort charge
#---------------
charge_all_new=[]
for i in range(0, len(Molecules)):
    charge=[]
    atom=atom_all[i]
    for j in atom:
        charge.append(str(float(j[6])*math.sqrt(138.935/er)))
    charge_all_new.append(charge)

#---------------
# type charge
#---------------
def exist1(array,name):
    for i in array:
        if i[0]==name:
            return True
    return False

charge_type=[]
for i in range(0, len(Molecules)):
    atom=atom_all[i]
    for j in range(0, len(atom)):
        if (not exist1(charge_type,atom[j][1])) and float(atom[j][6])!=0:
            charge_type.append([atom[j][1],float(atom[j][6])])



#---------------
# out put xml file
#---------------        
image=[]    
xml=open(grofile[0:len(grofile)-3]+"xml","w")   
print("output to "+grofile[0:len(grofile)-3]+"xml") 
xml.write('<?xml version ="1.0" encoding ="UTF-8" ?>\n')
xml.write('<galamost_xml version="1.3">\n')
xml.write('<configuration time_step="0" dimensions="3" natoms="'+str(len(pos))+'" >\n')
xml.write('<box lx="'+str(lx)+'" ly="'+str(ly)+'" lz="'+str(lz)+'"/>\n')
xml.write('<position num="'+str(len(pos))+'">\n')
for i in range(0, len(pos)):
    posix=pos[i][0]-lx/2.0
    posiy=pos[i][1]-ly/2.0
    posiz=pos[i][2]-lz/2.0
    ima=[0, 0, 0]
    if posix >lx/2.0:
        posix -= lx
        ima[0]=1
    elif posix <-lx/2.0:
        posix += lx
        ima[0]=-1
        
    if posiy >ly/2.0:
        posiy -= ly
        ima[1]=1
    elif posiy <-ly/2.0:
        posiy += ly
        ima[1]=-1

    if posiz >lz/2.0:
        posiz -= lz
        ima[2]=1
    elif posiz <-lz/2.0:
        posiz += lz
        ima[2]=-1
    image.append(ima)
    xml.write(str(posix)+"  "+str(posiy)+"  "+str(posiz)+"\n")
xml.write('</position>\n')

if len(vel) == len(pos):
    xml.write('<velocity num="'+str(len(vel))+'">\n')
    for i in range(0, len(vel)):
        xml.write(str(vel[i][0])+"  "+str(vel[i][1])+"  "+str(vel[i][2])+"\n")
    xml.write('</velocity>\n')

xml.write('<image num="'+str(len(image))+'">\n')
for i in range(0, len(image)):
    xml.write(str(image[i][0])+"  "+str(image[i][1])+"  "+str(image[i][2])+"\n")
xml.write('</image>\n') 

xml.write('<charge num="'+str(len(pos))+'">\n') 
count=0
for i in range(0, len(Molecules)):
    mol_num=int(Molecules[i][1])
    charge = charge_all_new[i]
    n=len(charge)
    for j in range(0, mol_num):
        for k in charge:
            xml.write(k+'\n')
        count += n
if len(pos) != count:
    print ("number of pos ",len(pos),", number of charge " ,count)
    raise RuntimeError('Error particle charge number is wrong!')
xml.write('</charge>\n')

xml.write('<type num="'+str(len(pos))+'">\n')   
count=0
for i in range(0, len(Molecules)):
    mol_num=int(Molecules[i][1])
    atom =atom_all_new[i]
    n=len(atom)
    for j in range(0, mol_num):
        for k in atom:
            xml.write(k+'\n')
        count += n
if len(pos) != count:
    raise RuntimeError('Error particle type number is wrong!')          
xml.write('</type>\n')

xml.write('<mass num="'+str(len(pos))+'">\n')   
count=0
for i in range(0, len(Molecules)):
    mol_num=int(Molecules[i][1])
    mass = mass_all_new[i]
    n=len(mass)
    for j in range(0, mol_num):
        for k in mass:
            xml.write(k+'\n')
        count += n
if len(pos) != count:
    raise RuntimeError('Error particle mass number is wrong!')          
xml.write('</mass>\n')

if len(bond_type)>0:
    xml.write('<bond>\n')
    count=0
    for i in range(0, len(Molecules)):
        mol_num=Molecules[i][1]
        atom =atom_all_new[i]
        n=len(atom)
        bond = bond_all_new[i]
        for j in range(0, mol_num):
            for k in bond:
                xml.write(k[0]+" "+str(k[1]+count)+" "+str(k[2]+count)+'\n')
            count+=n    
    xml.write('</bond>\n')

if len(constraint_type)>0:
    xml.write('<constraint>\n')
    count=0
    for i in range(0, len(Molecules)):
        mol_num=Molecules[i][1]
        atom =atom_all_new[i]
        n=len(atom)
        constraint = constraint_all_new[i]
        for j in range(0, mol_num):
            for k in constraint:
                xml.write(k[0]+" "+str(k[1]+count)+" "+str(k[2]+count)+'\n')
            count+=n    
    xml.write('</constraint>\n')

if len(virtual_type)>0:
    xml.write('<vsite>\n')
    count=0
    for i in range(0, len(Molecules)):
        mol_num=Molecules[i][1]
        atom =atom_all_new[i]
        n=len(atom) 
        virtual = virtual_all_new[i]
        for j in range(0, mol_num):
            for k in virtual:
                xml.write(k[0]+" "+str(k[1]+count)+" "+str(k[2]+count)+" "+str(k[3]+count)+" "+str(k[4]+count)+'\n')
            count+=n            
    xml.write('</vsite>\n')

if len(angle_type)>0:
    xml.write('<angle>\n')
    count=0
    for i in range(0, len(Molecules)):
        mol_num=Molecules[i][1]
        atom =atom_all_new[i]
        n=len(atom) 
        angle = angle_all_new[i]
        for j in range(0, mol_num):
            for k in angle:
                xml.write(k[0]+" "+str(k[1]+count)+" "+str(k[2]+count)+" "+str(k[3]+count)+'\n')
            count+=n            
    xml.write('</angle>\n')

if len(dihedral_type)>0:
    xml.write('<dihedral>\n')
    count=0
    for i in range(0, len(Molecules)):
        mol_num=Molecules[i][1]
        atom =atom_all_new[i]
        n=len(atom)     
        dihedral = dihedral_all_new[i]
        for j in range(0, mol_num):
            for k in dihedral:
                xml.write(k[0]+" "+str(k[1]+count)+" "+str(k[2]+count)+" "+str(k[3]+count)+" "+str(k[4]+count)+'\n')
            count+=n            
    xml.write('</dihedral>\n')

xml.write('</configuration>\n')
xml.write('</galamost_xml>\n')      
xml.close()

#---------------
# out put mol xml file
#---------------        
count=0
filelist=[]
for i in range(0, len(Molecules)):
    filename=Molecules[i][0]+".xml"
    mol_num=int(Molecules[i][1])
    atom =atom_all_new[i]
    np=len(atom)
    
    findname=False
    for j in filelist:
        if j == filename:
            count += mol_num*np
            findname=True
            break
            
    if findname:
        continue
    else:
        filelist.append(filename)
    
    xml = open(filename, "w")
    image=[]
    print(filename)   
    xml.write('<?xml version ="1.0" encoding ="UTF-8" ?>\n')
    xml.write('<galamost_xml version="1.3">\n')
    xml.write('<configuration time_step="0" dimensions="3" natoms="'+str(np)+'" >\n')
    xml.write('<box lx="'+str(lx)+'" ly="'+str(ly)+'" lz="'+str(lz)+'"/>\n')
    xml.write('<position num="'+str(np)+'">\n')
    for j in range(np):
        posix=pos[count+j][0]-lx/2.0
        posiy=pos[count+j][1]-ly/2.0
        posiz=pos[count+j][2]-lz/2.0
        ima=[0, 0, 0]
        if posix >lx/2.0:
            posix -= lx
            ima[0]=1
        elif posix <-lx/2.0:
            posix += lx
            ima[0]=-1
            
        if posiy >ly/2.0:
            posiy -= ly
            ima[1]=1
        elif posiy <-ly/2.0:
            posiy += ly
            ima[1]=-1
    
        if posiz >lz/2.0:
            posiz -= lz
            ima[2]=1
        elif posiz <-lz/2.0:
            posiz += lz
            ima[2]=-1
        image.append(ima)
        xml.write(str(posix)+"  "+str(posiy)+"  "+str(posiz)+"\n")
    xml.write('</position>\n')
    
    if len(vel) == len(pos):
        xml.write('<velocity num="'+str(np)+'">\n')
        for j in range(0, np):
            xml.write(str(vel[count+j][0])+"  "+str(vel[count+j][1])+"  "+str(vel[count+j][2])+"\n")
        xml.write('</velocity>\n')
    
    xml.write('<image num="'+str(len(image))+'">\n')
    for j in range(0, len(image)):
        xml.write(str(image[j][0])+"  "+str(image[j][1])+"  "+str(image[j][2])+"\n")
    xml.write('</image>\n') 
    
    xml.write('<charge num="'+str(np)+'">\n')
    charge = charge_all_new[i]
    for j in charge:
        xml.write(j+'\n')
    xml.write('</charge>\n')
    
    xml.write('<type num="'+str(np)+'">\n') 
    for j in atom:
        xml.write(j+'\n')
    xml.write('</type>\n')
    
    xml.write('<mass num="'+str(np)+'">\n')
    mass = mass_all_new[i]
    for j in mass:
        xml.write(j+'\n')
    xml.write('</mass>\n')
    
    bond = bond_all_new[i]
    if len(bond)>0:
        xml.write('<bond>\n')
        for j in bond:
            xml.write(j[0]+" "+str(j[1])+" "+str(j[2])+'\n')
        xml.write('</bond>\n')
    
    constraint = constraint_all_new[i]
    if len(constraint)>0:
        xml.write('<constraint>\n')
        for j in constraint:
            xml.write(j[0]+" "+str(j[1])+" "+str(j[2])+'\n')
        xml.write('</constraint>\n')
    
    virtual = virtual_all_new[i]
    if len(virtual)>0:
        xml.write('<vsite>\n')
        for j in virtual:
            xml.write(j[0]+" "+str(j[1])+" "+str(j[2])+" "+str(j[3])+" "+str(j[4])+'\n')
        xml.write('</vsite>\n')

    angle = angle_all_new[i]
    if len(angle)>0:
        xml.write('<angle>\n')
        for j in angle:
            xml.write(j[0]+" "+str(j[1])+" "+str(j[2])+" "+str(j[3])+'\n')
        xml.write('</angle>\n')
    
    dihedral = dihedral_all_new[i]
    if len(dihedral)>0:
        xml.write('<dihedral>\n')
        for j in dihedral:
            xml.write(j[0]+" "+str(j[1])+" "+str(j[2])+" "+str(j[3])+" "+str(j[4])+'\n')
        xml.write('</dihedral>\n')
        
    xml.write('</configuration>\n')
    xml.write('</galamost_xml>\n')
    xml.close()
    count += mol_num*np

###-------------parse non-bonded parameters
param=[]
value=[]
pair=[]
define=[]
defaults=[]


def existType(name):
    for i in atom_type:
        if i==name:
            return True
    return False

def ParseItp():
    read_defaults=False
    read_nonbond_params=False
    for filename in ItpFiles:
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
            if read_nonbond_params and len(la)<4:
                read_nonbond_params=False               
            if read_nonbond_params and len(la)==5 and existType(la[0]) and existType(la[1]):
                pair.append((la[0],la[1]))
                param.append((la[3],la[4]))
            if read_nonbond_params and len(la)==4 and existType(la[0]) and existType(la[1]):
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

            if exist(nospace,'[nonbond_params]') or exist(nospace,'[pairtypes]'):
                read_nonbond_params=True
#               print "found nonbond param"

        itp.close()

ParseItp()

if len(defaults)==0:
    defaults.append(("1","2"))

if len(param)!=0:
    for i in range(0,len(param)):
        c6=float(param[i][0])
        c12=float(param[i][1])
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
    
def get_sigma(name):
    for i in atomproperties:
        if i[0]==name:
            return float(i[len(i)-2])
def get_epsilon(name):
    for i in atomproperties:
        if i[0]==name:
            return float(i[len(i)-1])
        
if len(param)==0:
    for i in range(len(atom_type)):
        for j in range(i,len(atom_type)):
            pair.append([atom_type[i],atom_type[j]])
            if defaults[0][1]=='2':
                sig=(get_sigma(atom_type[i])+get_sigma(atom_type[j]))/2
                eps=(get_epsilon(atom_type[i])*get_epsilon(atom_type[j]))**0.5
                value.append([eps,sig])
            if defaults[0][1]=='3':
                sig=(get_sigma(atom_type[i])*get_sigma(atom_type[j]))**0.5
                eps=(get_epsilon(atom_type[i])*get_epsilon(atom_type[j]))**0.5
                value.append([eps,sig])

for i in range(0,len(atom_type)):
    for j in range(i,len(atom_type)):
        na0=atom_type[i]
        na1=atom_type[j]
        existed=False
        for p in pair:
            if na0==p[0] and na1==p[1]:
                existed=True
            elif na0==p[1] and na1==p[0]:
                existed=True
        if existed==False:
            print ("no exist",na0,na1)


ff=open(grofile[0:len(grofile)-3]+"force_field","w")    
print("output to "+grofile[0:len(grofile)-3]+"force_field") 


if len(atom_type)>0:
    ff.write('<pair_params>\n')
    for i in range(0,len(value)):
        outstr=pair[i][0]+" "+pair[i][1]+" "+str(value[i][0])+" "+str(value[i][1])+" 1.0"
        ff.write(outstr+'\n')
    ff.write('</pair_params>\n')
    ff.write('\n')  

if len(constraint_type)>0:
    ff.write('<constraint_params>\n')
    for i in range(0,len(constraint_type)):
        ff.write(constraint_type[i][0]+" "+constraint_type[i][2]+" "+constraint_type[i][1]+'\n')
    ff.write('</constraint_params>\n')
    ff.write('\n')  
        
if len(virtual_type)>0:
    ff.write('<vsite_params>\n')
    for i in range(0,len(virtual_type)):
        ff.write(virtual_type[i][0]+" "+virtual_type[i][2]+" "+virtual_type[i][3]+" "+virtual_type[i][4]+" "+virtual_type[i][1]+'\n')
    ff.write('</vsite_params>\n')   
    ff.write('\n')  
    
if len(bond_type)>0:
    ff.write('<bond_params>\n')
    for i in range(0,len(bond_type)):
        ff.write(bond_type[i][0]+" "+bond_type[i][3]+" "+bond_type[i][2]+" "+bond_type[i][1]+'\n')
    ff.write('</bond_params>\n')
    ff.write('\n')      

if len(angle_type)>0:
    ff.write('<angle_params>\n')
    for i in range(0,len(angle_type)):
        if len(angle_type[i])==4:
            ff.write(angle_type[i][0]+" "+angle_type[i][3]+" "+angle_type[i][2]+" "+angle_type[i][1]+'\n')
        elif len(angle_type[i])==6:
            ff.write(angle_type[i][0]+" "+angle_type[i][3]+" "+angle_type[i][2]+" "+angle_type[i][5]+" "+angle_type[i][4]+" "+angle_type[i][1]+'\n')
    ff.write('</angle_params>\n')
    ff.write('\n')  
    
if len(dihedral_type)>0:
    ff.write('<dihedral_params>\n')
    for i in range(0,len(dihedral_type)):
        ff.write(str(dihedral_type[i][0])+" ")
        for j in range(2, len(dihedral_type[i])):
            ff.write(str(dihedral_type[i][j])+" ")
        ff.write(str(dihedral_type[i][1])+" ")          
        ff.write('\n')
    ff.write('</dihedral_params>\n')
    
ff.close()



#---------------
# out put exclusion file
#---------------        
count=0
for i in range(0, len(Molecules)):
    eli = exclusion_all[i]
    count += len(eli)
if count > 0:
    counta = 0
    el = open("exclusion_list.dat", "w")
    print("output to exclusion_list.dat")
    for i in range(0, len(Molecules)):
        eli = exclusion_all[i]
        mol_num=int(Molecules[i][1])
        atom =atom_all_new[i]
        np=len(atom) 
        # print(Molecules[i][0], Molecules[i][1], len(eli))
        for j in range(0, mol_num):
            for k in eli:
                el.write(str(k[0]+counta)+" "+str(k[1]+counta)+"\n")
            counta += np
    el.close()
print ("success!")
            