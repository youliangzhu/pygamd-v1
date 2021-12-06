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

from pygamd.forces import potentials
import pygamd.snapshots.box as box_func
import numpy as np
import numba as nb
from numba import cuda
import math 
import time

def bond_force(cu_func):
	if isinstance(cu_func, str):
		cu_func = potentials.bonded_library.cu_bond(cu_func)
	if not isinstance(cu_func, cuda.compiler.DeviceDispatcher):
		raise RuntimeError('Error bond_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:, :], float32, float32, float32, float32, float32, float32, int32[:], int32[:, :, :], float32[:, :], float32[:, :], float32)")
	def cu_bond_force(npa, pos, params, box0, box1, box2, box0_half, box1_half, box2_half, bond_size, bond_list, force, virial_potential, one_sixth):
		i = cuda.grid(1)
		if i < npa:
			pix = pos[i][0]
			piy = pos[i][1]
			piz = pos[i][2]
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]	
			for bi in range(nb.int32(0), bond_size[i]):
				j = bond_list[i][bi][0]
				tj = bond_list[i][bi][1]
				dp0 = pix - pos[j][0]
				dp1 = piy - pos[j][1]
				dp2 = piz - pos[j][2]

				if dp0>=box0_half:
					dp0 -= box0
				elif dp0<-box0_half:
					dp0 += box0
					
				if dp1>=box1_half:
					dp1 -= box1
				elif dp1<-box1_half:
					dp1 += box1
				
				if dp2>=box2_half:
					dp2 -= box2
				elif dp2<-box2_half:
					dp2 += box2	

				rsq = dp0*dp0 + dp1*dp1 + dp2*dp2
				pms = params[tj]
				fp = cuda.local.array(2, dtype = nb.float32)
				cu_func(rsq, pms, fp)
				fijx = fp[0]*dp0
				fijy = fp[0]*dp1
				fijz = fp[0]*dp2
				virial = one_sixth * rsq * fp[0]
				potential = fp[1] * nb.float32(0.5)
				result0 += fijx
				result1 += fijy
				result2 += fijz
				result3 += virial
				result4 += potential
			force[i][0] = result0 
			force[i][1] = result1 
			force[i][2] = result2 
			virial_potential[i][0] = result3 
			virial_potential[i][1] = result4
	return cu_bond_force

class bond:
	#define init function
	def __init__(self, info, func):
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.bond.nbtypes)]
		self.block_size = 256
		self.bond_force = bond_force(func)
		self.params_changed = True
	#calculate bond force
	def calculate(self, timestep):
		if self.params_changed:
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False

		# for i in range(0, self.info.npa):
			# if self.info.bond.bond_size[i] >0:
				# print(i, self.info.bond.bond_size[i], self.info.bond.bond_list[i])
		# start = time.time()	
		self.bond_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
														self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.bond.d_bond_size,  
														self.info.bond.d_bond_list, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/6.0))
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)											   
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])

	def setParams(self, bond_type, param, check):
		type_id = self.info.bond.convert_name_to_id(bond_type)
		if isinstance(self.func, str) and self.func == 'harmonic':
			if len(param) != 2:
				raise RuntimeError('Error, the number of harmonic parameters is not equal to 2!')
			k=param[0]
			r0=param[1]
			self.params[type_id] = [k, r0]		
		elif isinstance(self.func, cuda.compiler.DeviceFunctionTemplate):
			self.params[type_id] = param

		check[type_id] = True
		self.params_changed = True
