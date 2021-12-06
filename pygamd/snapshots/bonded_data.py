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

import numpy as np
import numba as nb
from numba import cuda
import math

#sort bonds
@cuda.jit("void(int32, int32[:], int32[:, :, :], int32[:], int32[:, :, :], int32[:])")
def cu_sort_bond(npa, bond_size, bond_list, bond_size_tag, bond_list_tag, rtag):
	tag = cuda.grid(1)
	if tag < npa:
		nb = bond_size_tag[tag]
		idx = rtag[tag]
		for j in range(0, nb):
			bj = bond_list_tag[tag][j][0]
			tj = bond_list_tag[tag][j][1]
			
			bond_list[idx][j][0] = rtag[bj]
			bond_list[idx][j][1] = tj
		bond_size[idx] = nb


class bond_data:
	def __init__(self, info, bonds):
		self.info = info
		self.bonds = bonds
		self.bond_size = np.zeros(self.info.npa, dtype=np.int32)
		self.typemap=[]
		for i in range(0, len(bonds)):
			b = bonds[i]
			if b[1] >= self.info.npa or b[2] >= self.info.npa:
				raise RuntimeError('Error! particle idx ', b[1], b[2], ' in bonds exceed the number of particles ', self.npa)
			self.bond_size[b[1]] += np.int32(1)
			self.bond_size[b[2]] += np.int32(1)
			
		self.nmax_bond = np.int32(0)
		for i in range(0, self.info.npa):
			if self.bond_size[i] > self.nmax_bond:
				self.nmax_bond = self.bond_size[i]
			self.bond_size[i] = np.int32(0)
				
		self.bond_list = np.zeros([self.info.npa, self.nmax_bond, 2], dtype=np.int32)
		
		for i in range(0, len(bonds)):
			b = bonds[i]
			ti = self.add_name_to_id(b[0])
			
			nb1 = self.bond_size[b[1]]
			nb2 = self.bond_size[b[2]]
			
			self.bond_list[b[1]][nb1][0] = b[2]
			self.bond_list[b[2]][nb2][0] = b[1]
			
			self.bond_list[b[1]][nb1][1] = ti
			self.bond_list[b[2]][nb2][1] = ti		
			
			self.bond_size[b[1]] = nb1 + np.int32(1)
			self.bond_size[b[2]] = nb2 + np.int32(1)
			# print(b)
		# device arrays	
		self.d_bond_size = cuda.to_device(self.bond_size)
		self.d_bond_size_tag = cuda.to_device(self.bond_size)		
		self.d_bond_list = cuda.to_device(self.bond_list)		
		self.d_bond_list_tag = cuda.to_device(self.bond_list)
		
		self.nbtypes = len(self.typemap)
		self.block_size = 256
		self.nblocks = math.ceil(self.info.npa / self.block_size)
		
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
	# convert type name to type index		
	def sort(self):
		cu_sort_bond[self.nblocks, self.block_size](self.info.npa, self.d_bond_size, self.d_bond_list, self.d_bond_size_tag, self.d_bond_list_tag, self.info.d_rtag)


#sort angles
@cuda.jit("void(int32, int32[:], int32[:, :, :], int32[:], int32[:, :, :], int32[:])")
def cu_sort_angle(npa, angle_size, angle_list, angle_size_tag, angle_list_tag, rtag):
	tag = cuda.grid(1)
	if tag < npa:
		na = angle_size_tag[tag]
		idx = rtag[tag]
		for j in range(0, na):
			aj1 = angle_list_tag[tag][j][0]
			aj2 = angle_list_tag[tag][j][1]
			tj  = angle_list_tag[tag][j][2]
			order = angle_list_tag[tag][j][3]
			
			angle_list[idx][j][0] = rtag[aj1]
			angle_list[idx][j][1] = rtag[aj2]
			angle_list[idx][j][2] = tj
			angle_list[idx][j][3] = order			
		angle_size[idx] = na

class angle_data:
	def __init__(self, info, angles):
		self.info = info
		self.angles = angles
		self.angle_size = np.zeros(self.info.npa, dtype=np.int32)
		self.typemap=[]
		for i in range(0, len(angles)):
			a = angles[i]
			if a[1] >= self.info.npa or a[2] >= self.info.npa or a[3] >= self.info.npa:
				raise RuntimeError('Error! particle idx ', a[1], a[2], a[3], ' in angle exceed the number of particles ', self.npa)
			self.angle_size[a[1]] += np.int32(1)
			self.angle_size[a[2]] += np.int32(1)
			self.angle_size[a[3]] += np.int32(1)
			
		self.nmax_angle = np.int32(0)
		for i in range(0, self.info.npa):
			if self.angle_size[i] > self.nmax_angle:
				self.nmax_angle = self.angle_size[i]
			self.angle_size[i] = np.int32(0)
				
		self.angle_list = np.zeros([self.info.npa, self.nmax_angle, 4], dtype=np.int32)
		
		for i in range(0, len(angles)):
			a = angles[i]
			ti = self.add_name_to_id(a[0])
			
			na1 = self.angle_size[a[1]]
			na2 = self.angle_size[a[2]]
			na3 = self.angle_size[a[3]]			
			
			self.angle_list[a[1]][na1][0] = a[2]
			self.angle_list[a[2]][na2][0] = a[1]
			self.angle_list[a[3]][na3][0] = a[1]			
			
			self.angle_list[a[1]][na1][1] = a[3]
			self.angle_list[a[2]][na2][1] = a[3]
			self.angle_list[a[3]][na3][1] = a[2]

			self.angle_list[a[1]][na1][2] = ti
			self.angle_list[a[2]][na2][2] = ti
			self.angle_list[a[3]][na3][2] = ti

			self.angle_list[a[1]][na1][3] = np.int32(0)
			self.angle_list[a[2]][na2][3] = np.int32(1)
			self.angle_list[a[3]][na3][3] = np.int32(2)		
			
			self.angle_size[a[1]] = na1 + np.int32(1)
			self.angle_size[a[2]] = na2 + np.int32(1)
			self.angle_size[a[3]] = na3 + np.int32(1)			
			# print(a)
		# device arrays	
		self.d_angle_size = cuda.to_device(self.angle_size)
		self.d_angle_size_tag = cuda.to_device(self.angle_size)		
		self.d_angle_list = cuda.to_device(self.angle_list)		
		self.d_angle_list_tag = cuda.to_device(self.angle_list)
		
		self.natypes = len(self.typemap)
		self.block_size = 256
		self.nblocks = math.ceil(self.info.npa / self.block_size)
		
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
	# convert type name to type index		
	def sort(self):
		cu_sort_angle[self.nblocks, self.block_size](self.info.npa, self.d_angle_size, self.d_angle_list, self.d_angle_size_tag, self.d_angle_list_tag, self.info.d_rtag)



#sort dihedrals
@cuda.jit("void(int32, int32[:], int32[:, :, :], int32[:], int32[:, :, :], int32[:])")
def cu_sort_dihedral(npa, dihedral_size, dihedral_list, dihedral_size_tag, dihedral_list_tag, rtag):
	tag = cuda.grid(1)
	if tag < npa:
		nd = dihedral_size_tag[tag]
		idx = rtag[tag]
		for j in range(0, nd):
			aj1 = dihedral_list_tag[tag][j][0]
			aj2 = dihedral_list_tag[tag][j][1]
			aj3 = dihedral_list_tag[tag][j][2]			
			tj  = dihedral_list_tag[tag][j][3]
			order = dihedral_list_tag[tag][j][4]
			
			dihedral_list[idx][j][0] = rtag[aj1]
			dihedral_list[idx][j][1] = rtag[aj2]
			dihedral_list[idx][j][2] = rtag[aj3]			
			dihedral_list[idx][j][3] = tj
			dihedral_list[idx][j][4] = order			
		dihedral_size[idx] = nd

class dihedral_data:
	def __init__(self, info, dihedrals):
		self.info = info
		self.dihedrals = dihedrals
		self.dihedral_size = np.zeros(self.info.npa, dtype=np.int32)
		self.typemap=[]
		for i in range(0, len(dihedrals)):
			d = dihedrals[i]
			if d[1] >= self.info.npa or d[2] >= self.info.npa or d[3] >= self.info.npa or d[4] >= self.info.npa:
				raise RuntimeError('Error! particle idx ', d[1], d[2], d[3], d[4], ' in dihedral exceed the number of particles ', self.npa)
			self.dihedral_size[d[1]] += np.int32(1)
			self.dihedral_size[d[2]] += np.int32(1)
			self.dihedral_size[d[3]] += np.int32(1)
			self.dihedral_size[d[4]] += np.int32(1)
			
		self.nmax_dihedral = np.int32(0)
		for i in range(0, self.info.npa):
			if self.dihedral_size[i] > self.nmax_dihedral:
				self.nmax_dihedral = self.dihedral_size[i]
			self.dihedral_size[i] = np.int32(0)
				
		self.dihedral_list = np.zeros([self.info.npa, self.nmax_dihedral, 5], dtype=np.int32)
		
		for i in range(0, len(dihedrals)):
			d = dihedrals[i]
			ti = self.add_name_to_id(d[0])
			
			na1 = self.dihedral_size[d[1]]
			na2 = self.dihedral_size[d[2]]
			na3 = self.dihedral_size[d[3]]
			na4 = self.dihedral_size[d[4]]			
			
			self.dihedral_list[d[1]][na1][0] = d[2]
			self.dihedral_list[d[2]][na2][0] = d[1]
			self.dihedral_list[d[3]][na3][0] = d[1]
			self.dihedral_list[d[4]][na4][0] = d[1]				
			
			self.dihedral_list[d[1]][na1][1] = d[3]
			self.dihedral_list[d[2]][na2][1] = d[3]
			self.dihedral_list[d[3]][na3][1] = d[2]
			self.dihedral_list[d[4]][na4][1] = d[2]	

			self.dihedral_list[d[1]][na1][2] = d[4]
			self.dihedral_list[d[2]][na2][2] = d[4]
			self.dihedral_list[d[3]][na3][2] = d[4]
			self.dihedral_list[d[4]][na4][2] = d[3]				

			self.dihedral_list[d[1]][na1][3] = ti
			self.dihedral_list[d[2]][na2][3] = ti
			self.dihedral_list[d[3]][na3][3] = ti
			self.dihedral_list[d[4]][na4][3] = ti
			
			self.dihedral_list[d[1]][na1][4] = np.int32(0)
			self.dihedral_list[d[2]][na2][4] = np.int32(1)
			self.dihedral_list[d[3]][na3][4] = np.int32(2)		
			self.dihedral_list[d[4]][na4][4] = np.int32(3)
			
			self.dihedral_size[d[1]] = na1 + np.int32(1)
			self.dihedral_size[d[2]] = na2 + np.int32(1)
			self.dihedral_size[d[3]] = na3 + np.int32(1)
			self.dihedral_size[d[4]] = na4 + np.int32(1)			
			# print(d)
		# device arrays	
		self.d_dihedral_size = cuda.to_device(self.dihedral_size)
		self.d_dihedral_size_tag = cuda.to_device(self.dihedral_size)		
		self.d_dihedral_list = cuda.to_device(self.dihedral_list)		
		self.d_dihedral_list_tag = cuda.to_device(self.dihedral_list)
		
		self.ndtypes = len(self.typemap)
		self.block_size = 256
		self.nblocks = math.ceil(self.info.npa / self.block_size)
		
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
	# convert type name to type index		
	def sort(self):
		cu_sort_dihedral[self.nblocks, self.block_size](self.info.npa, self.d_dihedral_size, self.d_dihedral_list, self.d_dihedral_size_tag, self.d_dihedral_list_tag, self.info.d_rtag)


				