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

from pygamd import snapshots
import numpy as np
import numba as nb
from numba import cuda

# read a snapshot file with certain format
class read:
	def __init__(self, filename):
		file_array=filename.split(".")
		if file_array[len(file_array)-1]=="mst":
			self.data = snapshots.read_mst.read_mst(filename)
		self.compute_properties = {'temperature':False, 'pressure':False, 'momentum':False, 'potential':False, 'stress_tensor':False}
		self.variant = {'position':True, 'velocity':True, 'type':False, 'mass':False, 'image':True, 'box':False, 
						'force':True, 'potential':True, 'virial':True, 'bond':False, 'angle':False, 'dihedral':False}		
		# system information
		self.npa = self.data.num_particles
		self.pitch = (self.npa + (16 - (self.npa & 15)))
		self.dt=0.001
		self.timestep = self.data.timestep
		self.dimension = self.data.dimension
		
		self.particle_set = []
		self.comp_info = []
		self.plist = []
		self.typemap=[]
		
		self.charge = None
		self.body = None
		self.diameter = None
		self.molecule = None
		self.init = None		
		self.cris = None
		self.orientation = None
		self.quaternion = None	
		self.inert = None		
		
		
		self.bond = None
		self.angle = None
		self.dihedral = None		
		self.vsite = None
		
		# host arrays		
		self.pos = np.zeros([self.npa, 4], dtype=np.float32)
		self.vel = np.zeros([self.npa, 4], dtype=np.float32)
		self.image = np.zeros([self.npa, 3], dtype=np.int32)
		self.box = np.asarray(self.data.box, dtype=np.float32)
		self.tag = np.zeros(self.npa, dtype=np.int32)		
		self.rtag = np.zeros(self.npa, dtype=np.int32)		
		
		
		for i in range(0, self.npa):
			self.pos[i][0] = self.data.position[i][0]
			self.pos[i][1] = self.data.position[i][1]
			self.pos[i][2] = self.data.position[i][2]
			type_i = self.data.type[i]
			type_id = self.add_name_to_id(type_i)
			self.pos[i][3] = np.float32(type_id)
			self.tag[i] = i 
			self.rtag[i] = i

			
		for i in range(0, self.npa):
			if len(self.data.velocity)==self.npa:
				self.vel[i][0] = self.data.velocity[i][0]
				self.vel[i][1] = self.data.velocity[i][1]
				self.vel[i][2] = self.data.velocity[i][2]
			if len(self.data.mass)==self.npa:
				self.vel[i][3] = self.data.mass[i]
			else:
				self.vel[i][3] = 1.0		

		if len(self.data.image)==self.npa:
			for i in range(0, self.npa):
				self.image[i][0] = self.data.image[i][0]
				self.image[i][1] = self.data.image[i][1]
				self.image[i][2] = self.data.image[i][2]
				
		if len(self.data.charge)==self.npa:
			self.charge = np.zeros(self.npa, dtype=np.float32)
			for i in range(0, self.npa):
				self.charge[i] = self.data.charge[i]

		if len(self.data.body)==self.npa:
			self.body = np.zeros(self.npa, dtype=np.int32)
			for i in range(0, self.npa):
				self.body[i] = self.data.body[i]

		if len(self.data.diameter)==self.npa:
			self.diameter = np.zeros(self.npa, dtype=np.float32)
			for i in range(0, self.npa):
				self.diameter[i] = self.data.diameter[i]

		if len(self.data.molecule)==self.npa:
			self.molecule = np.zeros(self.npa, dtype=np.int32)
			for i in range(0, self.npa):
				self.molecule[i] = self.data.molecule[i]

		if len(self.data.init)==self.npa:
			self.init = np.zeros(self.npa, dtype=np.int32)
			for i in range(0, self.npa):
				self.init[i] = self.data.init[i]				

		if len(self.data.cris)==self.npa:
			self.cris = np.zeros(self.npa, dtype=np.int32)
			for i in range(0, self.npa):
				self.cris[i] = self.data.cris[i]
				
		if len(self.data.orientation)==self.npa:
			self.orientation = np.zeros([self.npa, 3], dtype=np.float32)
			for i in range(0, self.npa):
				self.orientation[i][0] = self.data.orientation[i][0]
				self.orientation[i][1] = self.data.orientation[i][1]
				self.orientation[i][2] = self.data.orientation[i][2]

		if len(self.data.quaternion)==self.npa:
			self.quaternion = np.zeros([self.npa, 4], dtype=np.float32)
			for i in range(0, self.npa):
				self.quaternion[i][0] = self.data.quaternion[i][0]
				self.quaternion[i][1] = self.data.quaternion[i][1]
				self.quaternion[i][2] = self.data.quaternion[i][2]
				self.quaternion[i][3] = self.data.quaternion[i][3]

		if len(self.data.inert)==self.npa:
			self.inert = np.zeros([self.npa, 3], dtype=np.float32)
			for i in range(0, self.npa):
				self.inert[i][0] = self.data.inert[i][0]
				self.inert[i][1] = self.data.inert[i][1]
				self.inert[i][2] = self.data.inert[i][2]
				
		self.force = np.zeros([self.npa, 3], dtype=np.float32)
		self.virial_potential = np.zeros([self.npa, 2], dtype=np.float32)
		self.ntypes=len(self.typemap)
		
		if len(self.data.bond)>0:
			self.bond = snapshots.bonded_data.bond_data(self, self.data.bond)
			
		if len(self.data.angle)>0:
			self.angle = snapshots.bonded_data.angle_data(self, self.data.angle)

		if len(self.data.dihedral)>0:
			self.dihedral = snapshots.bonded_data.dihedral_data(self, self.data.dihedral)			

		# device arrays		
		self.d_pos = cuda.to_device(self.pos)
		self.d_box = cuda.to_device(self.box)
		self.d_image = cuda.to_device(self.image)
		self.d_vel = cuda.to_device(self.vel)
		self.d_force = cuda.to_device(self.force)
		self.d_virial_potential = cuda.to_device(self.virial_potential)
		self.d_tag = cuda.to_device(self.tag)
		self.d_rtag = cuda.to_device(self.rtag)
		self.d_velo = None

		self.d_charge = None
		self.d_body = None
		self.d_diameter = None
		self.d_molecule = None
		self.d_init = None		
		self.d_cris = None
		self.d_orientation = None
		self.d_quaternion = None	
		self.d_inert = None

		if len(self.data.charge)==self.npa:
			self.d_charge = cuda.to_device(self.charge)

		if len(self.data.body)==self.npa:
			self.d_body = cuda.to_device(self.body)
			
		if len(self.data.diameter)==self.npa:
			self.d_diameter = cuda.to_device(self.diameter)

		if len(self.data.molecule)==self.npa:
			self.d_molecule = cuda.to_device(self.molecule)

		if len(self.data.init)==self.npa:
			self.d_init = cuda.to_device(self.init)			

		if len(self.data.cris)==self.npa:
			self.d_cris = cuda.to_device(self.cris)	
				
		if len(self.data.orientation)==self.npa:
			self.d_orientation = cuda.to_device(self.orientation)	

		if len(self.data.quaternion)==self.npa:
			self.d_quaternion = cuda.to_device(self.quaternion)	

		if len(self.data.inert)==self.npa:
			self.d_inert = cuda.to_device(self.inert)			

	# add type name and convert it to type index		
	def add_name_to_id(self, typename):
		for i in range(0, len(self.typemap)):
			if self.typemap[i] == typename:
				return i
		self.typemap.append(typename)
		return len(self.typemap)-1

	# convert type name to type index		
	def convert_name_to_id(self, typename):
		for i in range(0, len(self.typemap)):
			if self.typemap[i] == typename:
				return i
		raise RuntimeError('Error! type '+typename+' is not existed')
		

	def find_particle_set(self, group):
		for i in self.particle_set:
			if i.group==group:
				return i
			
	# find an existed comp_info with same group
	def find_comp_info(self, group):
		for i in self.comp_info:
			if i.group==group:
				return i
			
	def find_plist(self, rcut, exclusion):
		for i in self.plist:
			if i.rcut>=rcut and i.exclusion==exclusion:
				return i			
			
