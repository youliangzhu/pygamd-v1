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

from pygamd import chares
import numpy as np
from numba import cuda
import numba as nb
import math
		
class particle_set:
	def __init__(self, info, group):		
		self.info = info
		self.group = group
		self.check_member_tag = np.zeros(self.info.npa, dtype = np.int32)
		self.check_member_idx = np.zeros(self.info.npa, dtype = np.int32)		
		self.member=[]
		
		if isinstance(group, str):
			if group == "all":
				for i in range(0, self.info.npa):
					self.check_member_tag[i] = np.int32(1)
			else:
				raise RuntimeError("Error group key word '"+group+"'!")
		elif isinstance(group, list):
			for j in range(0, len(group)):
				element = group[j]
				if isinstance(element, str):
					count = 0
					for i in range(0, self.info.npa):
						typi = self.info.data.type[i]
						if typi == element:
							self.check_member_tag[i] = np.int32(1)
							count += 1
					if count ==0:
						 print ("Warnning, no particles with type '"+ element+"'")
				elif isinstance(element, int):
					self.check_member_tag[element] = np.int32(1)
				else:
					raise RuntimeError("Error type of group element in list ", type(element), "the type of elements should be string or int!")
		else:
			raise RuntimeError("Error group type ", type(group),", the type of group should be string or list !")
				
		for i in range(0, self.info.npa):
			if self.check_member_tag[i] == np.int32(1):
				self.member.append(i)	
	
		self.nme = len(self.member)
		if self.info.dimension == 3:
			self.nfreedom = self.nme*3-3			
		elif self.info.dimension == 2:
			self.nfreedom = self.nme*2-2		
		self.h_member = np.asarray(self.member, dtype=np.int32)
		self.d_member = cuda.to_device(self.h_member)
		
	def sort(self):
		self.info.tag = self.info.d_tag.copy_to_host()
		for i in range(0, self.info.npa):
			ti = self.info.tag[i]
			self.check_member_idx[i] = self.check_member_tag[ti]

		count = np.int32(0)
		for i in range(0, self.info.npa):
			if self.check_member_idx[i] == np.int32(1):
				self.h_member[count] = i
				count += np.int32(1)
		self.d_member = cuda.to_device(self.h_member)

		
class comp_info:
	def __init__(self, info, group):
		self.info = info
		self.group = group
		# self.ps = particle_set(info, group)
		
		self.ps = info.find_particle_set(group)
		if self.ps is None:
			self.ps = particle_set(info, group)
			info.particle_set.append(self.ps)		

		self.last_ts = 0xffffffff
		self.temp = 0.0
		self.pressure = 0.0	
		self.potential = 0.0
		self.momentum = 0.0
		self.block_size = 256
		self.nblocks = math.ceil(self.ps.nme / self.block_size)
		
		self.coll = np.zeros(self.nblocks*6, dtype = np.float32)
		self.d_coll = cuda.to_device(self.coll)
		# self.result = np.zeros(3, dtype = nb.float32)
		self.result = cuda.pinned_array(16, dtype = np.float32)
		self.result[0] = 0.0
		self.result[1] = 0.0
		self.result[2] = 0.0
		self.result[3] = 0.0
		self.result[4] = 0.0
		self.result[5] = 0.0		
		self.d_result = cuda.to_device(self.result)
		
	def calculate(self, timestep):
		if timestep == self.last_ts:
			return
		self.last_ts = timestep
		cp = self.info.compute_properties
		# print(timestep, cp)
		if cp['momentum'] and  cp['temperature'] and cp['pressure'] and cp['potential']:
			chares.compute_properties.cu_sums4[self.nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_vel, self.info.d_virial_potential, self.d_coll, self.nblocks)
			# h_coll=self.d_coll.copy_to_host()
			# print(self.nblocks)
			# print(h_coll)		
			chares.compute_properties.cu_final_sums4[1, 512](self.d_coll, self.d_result, self.nblocks)
		elif cp['temperature'] and cp['pressure']:
			chares.compute_properties.cu_sums2[self.nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_vel, self.info.d_virial_potential, self.d_coll, self.nblocks)
			chares.compute_properties.cu_final_sums2[1, 512](self.d_coll, self.d_result, self.nblocks)		
		elif cp['temperature']:
			chares.compute_properties.cu_sums1[self.nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_vel, self.info.d_virial_potential, self.d_coll, self.nblocks)
			chares.compute_properties.cu_final_sums1[1, 512](self.d_coll, self.d_result, self.nblocks)			
		
		
		self.result = self.d_result.copy_to_host()
		self.temp = self.result[0]/float(self.ps.nfreedom)
		
		volume_inv = 1.0/(self.info.box[0] * self.info.box[1] * self.info.box[2])
		virial = self.result[1]
		if self.info.dimension == 2:
			virial *= 3.0/2.0			
		self.pressure = (self.result[0] / float(self.info.dimension) + virial) * volume_inv
		self.potential = self.result[2]
		mx = self.result[3]
		my = self.result[4]
		mz = self.result[5]
		self.momentum = math.sqrt(mx*mx + my*my + mz*mz)/float(self.info.npa)
			# self.temp = sum_reduce(self.d_temp_arr)/float(self.nfreedom)
		# print (self.temp, self.result[1], self.result[2])
			
			
