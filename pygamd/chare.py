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

# [CuPy-native MVP] CuPy-native group membership + property reductions (no numba).
# particle_set keeps tag-defined membership on host for dump order, and maintains
# a CuPy current-index member list for reductions/integrators.  After particle
# sort, the device member list is rebuilt from d_tag on GPU (no tag host copy).
from pygamd import chares
from pygamd.backend.cupy import import_cupy
import numpy as np
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
		cp = import_cupy("CuPy particle_set")
		self.d_member = cp.asarray(self.h_member)
		self.d_check_member_tag = cp.asarray(self.check_member_tag)
		self._h_member_current = True

	def sort(self):
		# Rebuild the group's CURRENT-INDEX member list after a Sort reordered
		# particles.  Membership is defined by TAG (check_member_tag[tag]); after
		# reorder the particle at new index k has tag = d_tag[k], so it is a member
		# iff d_check_member_tag[d_tag[k]] == 1.  This stays on GPU.
		cp = import_cupy("CuPy particle_set")
		tag = self.info.cu.d_tag
		mask = self.d_check_member_tag[tag] != cp.int32(0)
		new_member = cp.nonzero(mask)[0].astype(cp.int32, copy=False)
		# Group size is tag-defined and should not change under sort.  Preserve the
		# d_member array identity when possible so graph signatures remain stable.
		if getattr(self, "d_member", None) is not None and self.d_member.shape == new_member.shape:
			self.d_member[...] = new_member
		else:
			self.d_member = new_member
		self.nme = int(self.d_member.shape[0])
		self._h_member_current = False

		
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
		# CuInfo CuPy device-data layer (required for the CuPy-native reductions).
		self.cu = getattr(info, "cu", None)
		if self.cu is None:
			raise NotImplementedError(
				"CuPy-native comp_info requires a CuInfo snapshot "
				"(pygamd.snapshot.read(..., backend='cupy')).")

	def calculate(self, timestep):
		if timestep == self.last_ts:
			return
		self.last_ts = timestep
		# CuPy-native group reduction (momentum / temperature / pressure / potential).
		ke2, virial, potential, mx, my, mz = chares.compute_properties.group_sums(
			self.cu, self.ps.d_member)
		self.temp = ke2 / float(self.ps.nfreedom)
		volume_inv = 1.0 / (self.info.box[0] * self.info.box[1] * self.info.box[2])
		if self.info.dimension == 2:
			virial *= 3.0 / 2.0
		self.pressure = (ke2 / float(self.info.dimension) + virial) * volume_inv
		self.potential = potential
		self.momentum = math.sqrt(mx * mx + my * my + mz * mz) / float(self.info.npa)
			
			
