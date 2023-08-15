'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with PYGAMD are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''
from numba.experimental import jitclass
from numba import types , typed
from numba.typed import Dict,List
from numba import char,int32,boolean,float32,float64
from numba import typeof
from numba import jit , objmode

#### ReadME
#### 1. Lines 470-560 have been commented as raise is not supported inside the objectmode. (https://numba.pydata.org/numba-doc/dev/user/withobjmode.html)
#### 2. 'asphere', 'patch','bond','angle','dihedral','vsite' have a string, try to convert it to integer like example conversion of 'type' done below.
#### 3. List of Lists are not supported, hence the code has been modified in snapshot.py as well.
#### 4. Below has implementations for Dict and List in Numba with jit support.


######### Numba throws errors with strings so this will pre-convert type to integer.
######### Function taken from snapshot.py
typemap = []

def add_name_to_id(typename):
		for i in range(0, len(typemap)):
			if typemap[i] == typename:
				return i
		typemap.append(typename)
		return len(typemap)-1
#########


######################### Numba supports only primitives types, but we can make it accept lists and Dictionaries using the following methods. 
######################### No guarantee with respect to performance though.

# Explicitly define the types of the key and value:
params_default = Dict.empty(key_type=typeof('box'),value_type=types.boolean)

# assign your default values
temp_dict = {'box':False, 'position':False, 'type':False, 'image':False, 'mass':False, 'velocity':False, 'charge':False, 'body':False,    
					   'diameter':False, 'rotangle':False, 'force':False, 'virial':False, 'molecule':False, 'init':False, 'cris':False, 
					   'orientation':False, 'quaternion':False, 'rotation':False, 'inert':False, 'asphere':False, 'patch':False,'bond':False, 
					   'angle':False, 'dihedral':False, 'vsite':False }	

# Host variable to copy inside the jit function
for x, y in temp_dict.items():
	params_default[x]=y

###### Explicit typing to let numba know what datatype the variables are. 
spec = [
    ('filename', char[:]),     # String
	('num_particles', int32),
	('timestep', int32),
	('line_no', int32),
	('dimension', boolean),
	('num_particles_read', boolean),
	('timestep_read', boolean),
	('dimension_read', boolean),
	('mst_read', boolean),
	('invariant_data', boolean),
	('variant_data', boolean),
	('read_indicator', typeof(params_default)),
]
	
# Creating List to store the values read from file. 	
for x, y in temp_dict.items():
	spec.append(tuple((str(x),types.ListType(float64))))
	spec.append(tuple((str(x)+"_read",boolean)))

# read a snapshot file with mst format
#### Using jitclass to accelerate all the underlying functions. ( https://numba.pydata.org/numba-doc/dev/user/jitclass.html )
@jitclass(spec)
class read_mst:
	# define init function
	def __init__(self,filename, params=params_default):
		self.num_particles=0
		self.timestep=0
		self.dimension=3
		self.line_no = 0
		self.mst_read =False
		self.invariant_data = False
		self.variant_data = False
		self.read_indicator = params
		self.init_data()	
		self.read_file(filename)
	
	def read_file(self, filename):
		with objmode(): ############# Important : objmode lets us run object code instead of compiled code in between jitted region. 
			file_object = open(filename)
			self.reset_params()
			for line in file_object:
				lin=line.strip('\n')
				line_array = lin.split()
				if len(line_array)==0:
					continue
					
				if line_array[0]=="mst_end":
					break

				if line_array[0]=="mst_version" and float(line_array[1])== 1.0:
					print("info : read mst file with version 1.0")
					self.mst_read = True
					continue
						
				if self.mst_read:
					if line_array[0]=="invariant_data":
						self.invariant_data = True
						self.variant_data = False						
						continue
						
					if line_array[0]=="variant_data":
						self.invariant_data = False
						self.variant_data = True	
						continue
						
					if line_array[0]=="frame":
						if self.variant_data:
							self.init_data()
						# else:
						# 	raise RuntimeError('Error! mst files with multiple frames without the label of "variant_data"')
						continue

					if line_array[0]=="num_particles":
						self.reset_params()
						self.num_particles_read=True
						continue
						
					if line_array[0]=="timestep":
						self.reset_params()
						self.timestep_read=True
						continue
						
					if line_array[0]=="dimension":
						self.reset_params()
						self.dimension_read=True
						continue						

					if line_array[0]=="box":
						self.reset_params()
						self.box_read=True
						if self.invariant_data:
							self.read_indicator['box'] = True
						continue

					if line_array[0]=="position":
						self.reset_params()				
						self.position_read=True
						if self.invariant_data:
							self.read_indicator['position'] = True						
						continue
					
					if line_array[0]=="type":
						self.reset_params()				
						self.type_read=True
						if self.invariant_data:
							self.read_indicator['type'] = True						
						continue
						
					if line_array[0]=="image":
						self.reset_params()
						self.image_read=True
						if self.invariant_data:
							self.read_indicator['image'] = True
						continue
						
					if line_array[0]=="mass":
						self.reset_params()
						self.mass_read=True
						if self.invariant_data:
							self.read_indicator['mass'] = True						
						continue
						
					if line_array[0]=="velocity":
						self.reset_params()
						self.velocity_read=True
						if self.invariant_data:
							self.read_indicator['velocity'] = True
						continue
						
					if line_array[0]=="charge":
						self.reset_params()
						self.charge_read=True
						if self.invariant_data:
							self.read_indicator['charge'] = True
						continue
						
					if line_array[0]=="body":
						self.reset_params()
						self.body_read=True
						if self.invariant_data:
							self.read_indicator['body'] = True
						continue

					if line_array[0]=="diameter":
						self.reset_params()
						self.diameter_read=True
						if self.invariant_data:
							self.read_indicator['diameter'] = True
						continue	

					if line_array[0]=="rotangle":
						self.reset_params()
						self.rotangle_read=True
						if self.invariant_data:
							self.read_indicator['rotangle'] = True
						continue
						
					if line_array[0]=="force":
						self.reset_params()
						self.force_read=True
						if self.invariant_data:
							self.read_indicator['force'] = True
						continue

					if line_array[0]=="virial":
						self.reset_params()
						self.virial_read=True
						if self.invariant_data:
							self.read_indicator['virial'] = True
						continue						

					if line_array[0]=="molecule":
						self.reset_params()
						self.molecule_read=True
						if self.invariant_data:
							self.read_indicator['molecule'] = True
						continue

					if line_array[0]=="init":
						self.reset_params()
						self.init_read=True
						if self.invariant_data:
							self.read_indicator['init'] = True
						continue	

					if line_array[0]=="cris":
						self.reset_params()
						self.cris_read=True
						if self.invariant_data:
							self.read_indicator['cris'] = True
						continue	

					if line_array[0]=="orientation":
						self.reset_params()
						self.orientation_read=True
						if self.invariant_data:
							self.read_indicator['orientation'] = True
						continue	

					if line_array[0]=="quaternion":
						self.reset_params()
						self.quaternion_read=True
						if self.invariant_data:
							self.read_indicator['quaternion'] = True
						continue	

					if line_array[0]=="rotation":
						self.reset_params()
						self.rotation_read=True
						if self.invariant_data:
							self.read_indicator['rotation'] = True
						continue	

					if line_array[0]=="inert":
						self.reset_params()
						self.inert_read=True
						if self.invariant_data:
							self.read_indicator['inert'] = True
						continue						

					if line_array[0]=="asphere":
						self.reset_params()
						self.asphere_read=True	
						if self.invariant_data:
							self.read_indicator['asphere'] = True
						continue	

					if line_array[0]=="patch":
						self.reset_params()
						self.patch_read=True	
						if self.invariant_data:
							self.read_indicator['patch'] = True
						continue		

					if line_array[0]=="bond":
						self.reset_params()
						self.bond_read=True	
						if self.invariant_data:
							self.read_indicator['bond'] = True
						continue						
						
					if line_array[0]=="angle":
						self.reset_params()
						self.angle_read=True
						if self.invariant_data:
							self.read_indicator['angle'] = True
						continue
						
					if line_array[0]=="dihedral":
						self.reset_params()
						self.dihedral_read=True
						if self.invariant_data:
							self.read_indicator['dihedral'] = True
						continue					
				
					if line_array[0]=="vsite":
						self.reset_params()
						self.vsite_read=True
						if self.invariant_data:
							self.read_indicator['vsite'] = True
						continue					
						
					# read data
					if self.num_particles_read and len(line_array) == 1:
						self.num_particles=int(line_array[0])
						print("info : number of particles", self.num_particles)
					
					if self.timestep_read and len(line_array) == 1:
						self.timestep=int(line_array[0])
						print("info : timestep", self.timestep)
						
					if self.dimension_read and len(line_array) == 1:
						self.dimension=int(line_array[0])		
						print("info : dimension", self.dimension)
						
					if self.box_read and len(line_array) == 3:
						self.box.append(float(line_array[0]))
						self.box.append(float(line_array[1]))
						self.box.append(float(line_array[2]))
						print("info : box size", line_array[0], line_array[1], line_array[2])
						
					if self.position_read and len(line_array) == 3:
						self.position.append(float(line_array[0]))
						self.position.append(float(line_array[1]))
						self.position.append(float(line_array[2]))
						
					if self.type_read and len(line_array) == 1:
						self.type.append(float(add_name_to_id(line_array[0])))
						
					if self.image_read and len(line_array) == 3:
						self.image.append(int(line_array[0]))
						self.image.append(int(line_array[1]))
						self.image.append(int(line_array[2]))
						
					if self.mass_read and len(line_array) == 1:
						self.mass.append(float(line_array[0]))
						
					if self.velocity_read and len(line_array) == 3:
						self.velocity.append(float(line_array[0]))
						self.velocity.append(float(line_array[1]))
						self.velocity.append(float(line_array[2]))
						
					if self.charge_read and len(line_array) == 1:
						self.charge.append(float(line_array[0]))

					if self.body_read and len(line_array) == 1:
						self.body.append(int(line_array[0]))

					if self.diameter_read and len(line_array) == 1:
						self.diameter.append(float(line_array[0]))

					if self.rotangle_read and len(line_array) == 3:
						self.rotangle.append(float(line_array[0]))
						self.rotangle.append(float(line_array[1]))
						self.rotangle.append(float(line_array[2]))							

					if self.force_read and len(line_array) == 3:
						self.force.append(float(line_array[0]))
						self.force.append(float(line_array[1]))
						self.force.append(float(line_array[2]))
						
					if self.virial_read and len(line_array) == 1:
						self.virial.append(float(line_array[0]))	

					if self.molecule_read and len(line_array) == 1:
						self.molecule.append(int(line_array[0]))

					if self.init_read and len(line_array) == 1:
						self.init.append(int(line_array[0]))

					if self.cris_read and len(line_array) == 1:
						self.cris.append(int(line_array[0]))

					if self.orientation_read and len(line_array) == 3:
						self.orientation.append(float(line_array[0]))
						self.orientation.append(float(line_array[1]))
						self.orientation.append(float(line_array[2]))						

					if self.quaternion_read and len(line_array) == 4:
						self.quaternion.append(float(line_array[0]))
						self.quaternion.append(float(line_array[1]))
						self.quaternion.append(float(line_array[2]))
						self.quaternion.append(float(line_array[3]))
						
					if self.rotation_read and len(line_array) == 3:
						self.rotation.append(float(line_array[0]))
						self.rotation.append(float(line_array[1]))
						self.rotation.append(float(line_array[2]))

					if self.inert_read and len(line_array) == 3:
						self.inert.append(float(line_array[0]))
						self.inert.append(float(line_array[1]))
						self.inert.append(float(line_array[2]))

					# if self.asphere_read and len(line_array) == 7:
					# 	self.asphere.append(line_array[0])
					# 	self.asphere.append(float(line_array[1]))
					# 	self.asphere.append(float(line_array[2]))
					# 	self.asphere.append(float(line_array[3]))
					# 	self.asphere.append(float(line_array[4]))
					# 	self.asphere.append(float(line_array[5]))
					# 	self.asphere.append(float(line_array[6]))						

					# if self.patch_read and len(line_array) == 5:
					# 	self.patch.append(line_array[0])
					# 	self.patch.append(float(line_array[1]))
					# 	self.patch.append(float(line_array[2]))
					# 	self.patch.append(float(line_array[3]))
					# 	self.patch.append(float(line_array[4]))						

					# if self.bond_read and len(line_array) == 3:
					# 	self.bond.append(line_array[0])
					# 	self.bond.append(int(line_array[1]))
					# 	self.bond.append(int(line_array[2]))

					# if self.angle_read and len(line_array) == 4:
					# 	self.angle.append(line_array[0])
					# 	self.angle.append(int(line_array[1]))
					# 	self.angle.append(int(line_array[2]))
					# 	self.angle.append(int(line_array[3]))

					# if self.dihedral_read and len(line_array) == 5:
					# 	self.dihedral.append(line_array[0])
					# 	self.dihedral.append(int(line_array[1]))
					# 	self.dihedral.append(int(line_array[2]))
					# 	self.dihedral.append(int(line_array[3]))
					# 	self.dihedral.append( int(line_array[4]))

					# if self.vsite_read and len(line_array) == 5:
					# 	self.vsite.append(line_array[0])
					# 	self.vsite.append(int(line_array[1]))
					# 	self.vsite.append(int(line_array[2]))
					# 	self.vsite.append(int(line_array[3]))
					# 	self.vsite.append(int(line_array[4]))
			file_object.close()						
						
		# # check data
		# if not self.mst_read:
		# 	raise RuntimeError('Error! the input file is not a mst file with version 1.0!')
		
		# if len(self.position) != self.num_particles:
		# 	raise RuntimeError('Error! number of position ', len(self.position), ' is not equal to the number of particles ', self.num_particles)	
		# print("info :", len(self.position), "positions")

		# if len(self.type) != self.num_particles:
		# 	raise RuntimeError('Error! number of type ', len(self.type), ' is not equal to the number of particles ', self.num_particles)
		# print("info :", len(self.type), "types")
		
		# if len(self.image) > 0:
		# 	if len(self.image) != self.num_particles:
		# 		raise RuntimeError('Error! number of image ', len(self.image), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.image), "images")
					
		# if len(self.mass) > 0:
		# 	if len(self.mass) != self.num_particles:
		# 		raise RuntimeError('Error! number of mass ', len(self.mass), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.mass), "masses")
			
		# if len(self.velocity) > 0:
		# 	if len(self.velocity) != self.num_particles:
		# 		raise RuntimeError('Error! number of velocity ', len(self.velocity), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.velocity), "velocities")
			
		# if len(self.charge) > 0:
		# 	if len(self.charge) != self.num_particles:
		# 		raise RuntimeError('Error! number of charge ', len(self.charge), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.charge), "charges")

		# if len(self.body) > 0:
		# 	if len(self.body) != self.num_particles:
		# 		raise RuntimeError('Error! number of body ', len(self.body), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.body), "bodies")

		# if len(self.diameter) > 0:
		# 	if len(self.diameter) != self.num_particles:
		# 		raise RuntimeError('Error! number of diameter ', len(self.diameter), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.diameter), "diameters")

		# if len(self.rotangle) > 0:
		# 	if len(self.rotangle) != self.num_particles:
		# 		raise RuntimeError('Error! number of rotangle ', len(self.rotangle), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.rotangle), "rotangles")			
			
		# if len(self.force) > 0:
		# 	if len(self.force) != self.num_particles:
		# 		raise RuntimeError('Error! number of force ', len(self.force), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.force), "forces")

		# if len(self.virial) > 0:
		# 	if len(self.virial) != self.num_particles:
		# 		raise RuntimeError('Error! number of virial ', len(self.virial), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.virial), "virials")

		# if len(self.molecule) > 0:
		# 	if len(self.molecule) != self.num_particles:
		# 		raise RuntimeError('Error! number of molecule ', len(self.molecule), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.molecule), "molecules")

		# if len(self.init) > 0:
		# 	if len(self.init) != self.num_particles:
		# 		raise RuntimeError('Error! number of init ', len(self.init), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.init), "inits")

		# if len(self.cris) > 0:
		# 	if len(self.cris) != self.num_particles:
		# 		raise RuntimeError('Error! number of cris ', len(self.cris), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.cris), "crises")

		# if len(self.orientation) > 0:
		# 	if len(self.orientation) != self.num_particles:
		# 		raise RuntimeError('Error! number of orientation ', len(self.orientation), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.orientation), "orientations")

		# if len(self.quaternion) > 0:
		# 	if len(self.quaternion) != self.num_particles:
		# 		raise RuntimeError('Error! number of quaternion ', len(self.quaternion), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.quaternion), "quaternions")

		# if len(self.rotation) > 0:
		# 	if len(self.rotation) != self.num_particles:
		# 		raise RuntimeError('Error! number of rotation ', len(self.rotation), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.rotation), "rotations")

		# if len(self.inert) > 0:
		# 	if len(self.inert) != self.num_particles:
		# 		raise RuntimeError('Error! number of inert ', len(self.inert), ' is not equal to the number of particles ', self.num_particles)
		# 	print("info :", len(self.inert), "inerts")			

		if len(self.asphere) > 0:
			print("info :", len(self.asphere), "aspheres")
			
		if len(self.patch) > 0:
			print("info :", len(self.patch), "patches")			
			
		if len(self.bond) > 0:
			print("info :", len(self.bond), "bonds")
			
		if len(self.angle) > 0:
			print("info :", len(self.angle), "angles")

		if len(self.dihedral) > 0:
			print("info :", len(self.dihedral), "dihedrals")

		if len(self.vsite) > 0:
			print("info :", len(self.vsite), "vsites")

	def init_data(self):
		# data
		if not self.read_indicator['box']:				
			self.box= typed.List.empty_list(types.float64)
		if not self.read_indicator['position']:				
			self.position= typed.List.empty_list(types.float64)			
		if not self.read_indicator['type']:			
			self.type= typed.List.empty_list(types.float64)
		if not self.read_indicator['image']:			
			self.image= typed.List.empty_list(types.float64)
		if not self.read_indicator['mass']:			
			self.mass= typed.List.empty_list(types.float64)
		if not self.read_indicator['velocity']:
			self.velocity= typed.List.empty_list(types.float64)
		if not self.read_indicator['charge']:
			self.charge = typed.List.empty_list(types.float64)
		if not self.read_indicator['body']:
			self.body= typed.List.empty_list(types.float64)
		if not self.read_indicator['diameter']:
			self.diameter= typed.List.empty_list(types.float64)
		if not self.read_indicator['rotangle']:
			self.rotangle= typed.List.empty_list(types.float64)				
		if not self.read_indicator['force']:			
			self.force= typed.List.empty_list(types.float64)
		if not self.read_indicator['virial']:
			self.virial= typed.List.empty_list(types.float64)
		if not self.read_indicator['molecule']:
			self.molecule= typed.List.empty_list(types.float64)
		if not self.read_indicator['init']:
			self.init= typed.List.empty_list(types.float64)
		if not self.read_indicator['cris']:
			self.cris= typed.List.empty_list(types.float64)
		if not self.read_indicator['orientation']:
			self.orientation= typed.List.empty_list(types.float64)
		if not self.read_indicator['quaternion']:
			self.quaternion= typed.List.empty_list(types.float64)
		if not self.read_indicator['rotation']:
			self.rotation= typed.List.empty_list(types.float64)			
		if not self.read_indicator['inert']:
			self.inert= typed.List.empty_list(types.float64)
		if not self.read_indicator['asphere']:
			self.asphere= typed.List.empty_list(types.float64)
		if not self.read_indicator['patch']:
			self.patch= typed.List.empty_list(types.float64)			
		if not self.read_indicator['bond']:
			self.bond= typed.List.empty_list(types.float64)			
		if not self.read_indicator['angle']:
			self.angle= typed.List.empty_list(types.float64)
		if not self.read_indicator['dihedral']:
			self.dihedral= typed.List.empty_list(types.float64)
		if not self.read_indicator['vsite']:
			self.vsite= typed.List.empty_list(types.float64)
			
	# reset parameters
	def reset_params(self):
		# indicators
		self.num_particles_read=False
		self.timestep_read=False
		self.dimension_read=False		
		self.box_read=False			
		self.position_read=False
		self.type_read=False
		self.image_read=False
		self.mass_read=False
		self.velocity_read=False
		self.charge_read=False	
		self.body_read=False
		self.diameter_read=False
		self.rotangle_read=False
		self.force_read=False
		self.virial_read=False
		self.molecule_read=False
		self.init_read=False
		self.cris_read=False
		self.orientation_read=False
		self.quaternion_read=False
		self.rotation_read=False
		self.inert_read=False
		self.asphere_read=False	
		self.patch_read=False			
		self.bond_read=False
		self.angle_read=False
		self.dihedral_read=False
		self.vsite_read=False