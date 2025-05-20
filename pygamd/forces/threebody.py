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

# try:
	# deviceFunction = cuda.compiler.DeviceFunctionTemplate
# except:
	# deviceFunction = cuda.compiler.DeviceDispatcher

# non-bonded force for neutral particles
def nonbonded_tri_force(cu_func):
	# if isinstance(cu_func, str):
		# cu_func = potentials.nb_library.cu_nonbonded(cu_func)
	# if not isinstance(cu_func, deviceFunction):
		# raise RuntimeError('Error nonbonded_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:, :], float32, float32, float32, int32, int32[:], int32[:, :], float32[:, :], float32[:, :])")
	def cu_nonbonded_tri_force(npa, pos, params, box0_half, box1_half, box2_half, ntypes, neighbor_size, neighbor_list, force, virial_potential):
		i = cuda.grid(1)
		if i < npa:
			pi = pos[i]
			ti = nb.int32(pi[3])
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]
			for ni in range(nb.int32(0), neighbor_size[i]):
				j = neighbor_list[i][ni]
				pj = pos[j]
				tj = nb.int32(pj[3])
				for nj in range(ni+1, neighbor_size[i]):
					k = neighbor_list[i][nj]
					pk = pos[k]
					tk = nb.int32(pk[3])

					pms = params[ti + tj*ntypes + tk*ntypes*ntypes]
					fp = cuda.local.array(5, dtype = nb.float32)
					cu_func(pi, pj, pk, box0_half, box1_half, box2_half, pms, fp)
					fijx = fp[0]
					fijy = fp[1]
					fijz = fp[2]
					virial = fp[3]
					potential = fp[4]
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
	return cu_nonbonded_tri_force

class threebody:
	#define init function
	def __init__(self, info, nlist, func):
		self.nlist=nlist
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.ntypes*self.info.ntypes*self.info.ntypes)]
		self.block_size = 256
		self.nonbonded_tri_force = nonbonded_tri_force(func)
		self.params_changed = True
	#calculate non-bonded force
	def calculate(self, timestep):
		if self.params_changed:
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False
			
		self.nlist.calculate(timestep)
		d_nlist = self.nlist.data.d_neighbor_list
		d_nlist_size = self.nlist.data.d_neighbor_size 
		
		# test
		# neigh_size = np.zeros([self.info.npa], dtype = np.int32)
		# neigh_size[1] = 5
		# neigh_list = np.zeros([self.info.npa, 8], dtype = np.int32)
		# neigh_list[1][0] = 0
		# neigh_list[1][1] = 4
		# neigh_list[1][2] = 5
		# neigh_list[1][3] = 2
		# neigh_list[1][4] = 6
		# d_nlist_size = cuda.to_device(neigh_size)
		# d_nlist = cuda.to_device(neigh_list)		
		# start = time.time()	
		self.nonbonded_tri_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_params,
															self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.ntypes, d_nlist_size, 
															d_nlist, self.info.d_force, self.info.d_virial_potential)
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
				
	#calculate non-bonded force in host
	# def calculate_host(self, timestep):
#		print(self.info.virial_potential)

	def setParams(self, type_i, type_j, type_k, param, check):
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)
		type_k_id = self.info.convert_name_to_id(type_k)
		idx1 = type_i_id * self.info.ntypes * self.info.ntypes + type_j_id * self.info.ntypes + type_k_id
		idx2 = type_k_id * self.info.ntypes * self.info.ntypes + type_j_id * self.info.ntypes + type_i_id
		self.params[idx1] = param
		self.params[idx2] = param
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True



