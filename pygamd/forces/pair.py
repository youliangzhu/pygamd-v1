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
def nonbonded_force(cu_func):
	if isinstance(cu_func, str):
		cu_func = potentials.nb_library.cu_nonbonded(cu_func)
	# if not isinstance(cu_func, deviceFunction):
		# raise RuntimeError('Error nonbonded_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:, :], float32, float32, float32, float32, float32, float32, int32, int32[:], int32[:, :], float32[:, :], float32[:, :], float32)")
	def cu_nonbonded_force(npa, pos, params, box0, box1, box2, box0_half, box1_half, box2_half, ntypes, neighbor_size, neighbor_list, force, virial_potential, one_sixth):
		i = cuda.grid(1)
		if i < npa:
			pix = pos[i][0]
			piy = pos[i][1]
			piz = pos[i][2]			
			ti = nb.int32(pos[i][3])
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]	
			for ni in range(nb.int32(0), neighbor_size[i]):
				j = neighbor_list[i][ni]
				# pj = pos[j]
				# dp = cuda.local.array(3, dtype = nb.float32)
				dp0 = pix - pos[j][0]
				dp1 = piy - pos[j][1]
				dp2 = piz - pos[j][2]
				tj = nb.int32(pos[j][3])
				# box_func.cu_box_min_dis(dp, box)
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
				pms = params[ti + tj*ntypes]
				fp = cuda.local.array(2, dtype = nb.float32)
				fp[0] = np.float32(0.0)
				fp[1] = np.float32(0.0)				
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
	return cu_nonbonded_force

class pair:
	#define init function
	def __init__(self, info, nlist, func):
		self.nlist=nlist
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.ntypes*self.info.ntypes)]
		self.block_size = 256
		if isinstance(self.func, str) and self.func == 'lj_coulomb':
			raise RuntimeError('Error, lj_coulomb should be called by nonbonded_c, not nonbonded!')
		self.nonbonded_force = nonbonded_force(func)
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
		self.nonbonded_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
															self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.ntypes, d_nlist_size, 
															d_nlist, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/6.0))
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)											   
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
				
	#calculate non-bonded force in host
	def calculate_host(self, timestep):
		npa = self.info.npa
		pos = self.info.pos
		typ = self.info.type
		self.nlist.calculate(timestep)
		list= self.nlist.data.list
		nb_func = potentials.nb_library.nonbonded(self.func)
		for idx in range(0, npa):
			nn=len(list[idx])
			pi = pos[idx]
			type_i = typ[idx]
			result=np.array([[0.0, 0.0, 0.0], [0.0, 0.0]])
			for j in range(0, nn):
				neighbor = list[idx][j]
				pj=pos[neighbor]
				type_j = typ[neighbor]
				
				dp = pi - pj
				box_func.box_min_dis(dp, self.info.box)
				param = self.params[type_i*self.info.ntypes + type_j]
				fp=[0.0, 0.0]
				rsq=(dp*dp).sum()
#				print(self.func)
				nb_func(rsq, param, fp)
				fij = fp[0]*dp
				result[0] += fij
				result[1][0] += (fij*dp).sum()
				result[1][1] += fp[1]
			forcei = self.info.force[idx]
			virial_potentiali = self.info.virial_potential[idx]
			self.info.force[idx] = forcei+result[0]
			self.info.virial_potential[idx] = virial_potentiali+result[1]
#		print(self.info.virial_potential)

	def setParams(self, type_i, type_j, param, check):	
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)	
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id
		# print(type_i_id, type_j_id, idx1, idx2, self.info.ntypes)
		if isinstance(self.func, str) and self.func == 'lj':
			if len(param) != 4:
				raise RuntimeError('Error, the number of lj parameters is not equal to 4!')
			epsilon=param[0]
			sigma=param[1]
			alpha=param[2]
			rcut=param[3]
			p0 = 4.0 * epsilon * math.pow(sigma, int(12))
			p1 = alpha * 4.0 * epsilon * math.pow(sigma, int(6))
			rcutsq = rcut*rcut
			self.params[idx1] = [p0, p1, rcutsq]
			self.params[idx2] = [p0, p1, rcutsq]
		elif isinstance(self.func, str) and self.func == 'harmonic':
			if len(param) != 2:
				raise RuntimeError('Error, the number of harmonic parameters is not equal to 2!')
			alpha=param[0]
			rcut=param[1]
			rcutsq = rcut*rcut
			self.params[idx1] = [alpha, 1.0/rcut, rcutsq]
			self.params[idx2] = [alpha, 1.0/rcut, rcutsq]
		else:
			self.params[idx1] = param
			self.params[idx2] = param
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True


# non-bonded force for charge particles
def nonbonded_c_force(cu_func):
	if isinstance(cu_func, str):
		cu_func = potentials.nb_library.cu_nonbonded(cu_func)
	# if not isinstance(cu_func, deviceFunction):
		# raise RuntimeError('Error nonbonded_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:], float32[:, :], float32, float32, float32, float32, float32, float32, int32, int32[:], int32[:, :], float32[:, :], float32[:, :], float32)")
	def cu_nonbonded_c_force(npa, pos, charge, params, box0, box1, box2, box0_half, box1_half, box2_half, ntypes, neighbor_size, neighbor_list, force, virial_potential, one_sixth):
		i = cuda.grid(1)
		if i < npa:
			pix = pos[i][0]
			piy = pos[i][1]
			piz = pos[i][2]
			qi = charge[i]
			ti = nb.int32(pos[i][3])
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]	
			for ni in range(nb.int32(0), neighbor_size[i]):
				j = neighbor_list[i][ni]
				# pj = pos[j]
				# dp = cuda.local.array(3, dtype = nb.float32)
				dp0 = pix - pos[j][0]
				dp1 = piy - pos[j][1]
				dp2 = piz - pos[j][2]
				qj = charge[j]
				tj = nb.int32(pos[j][3])
				# box_func.cu_box_min_dis(dp, box)
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
				pms = params[ti + tj*ntypes]
				fp = cuda.local.array(2, dtype = nb.float32)
				fp[0] = np.float32(0.0)
				fp[1] = np.float32(0.0)
				cu_func(rsq, qi, qj, pms, fp)
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
	return cu_nonbonded_c_force

class pair_c:
	#define init function
	def __init__(self, info, nlist, func):
		self.nlist=nlist
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.ntypes*self.info.ntypes)]
		self.block_size = 256
		if isinstance(self.func, str) and self.func == 'lj':
			raise RuntimeError('Error, lj should be called by nonbonded, not nonbonded_c!')
		elif isinstance(self.func, str) and self.func == 'harmonic':
			raise RuntimeError('Error, harmonic should be called by nonbonded, not nonbonded_c!')        
		self.nonbonded_c_force = nonbonded_c_force(func)
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
		self.nonbonded_c_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.info.d_charge, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
															self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.ntypes, d_nlist_size, 
															d_nlist, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/6.0))
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)											   
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
				
	#calculate non-bonded force in host
	def calculate_host(self, timestep):
		npa = self.info.npa
		pos = self.info.pos
		charge = self.info.charge
		typ = self.info.type
		self.nlist.calculate(timestep)
		list= self.nlist.data.list
		nb_func = potentials.nb_library.nonbonded(self.func)
		for idx in range(0, npa):
			nn=len(list[idx])
			pi = pos[idx]
			qi = charge[idx]
			type_i = typ[idx]
			result=np.array([[0.0, 0.0, 0.0], [0.0, 0.0]])
			for j in range(0, nn):
				neighbor = list[idx][j]
				pj = pos[neighbor]
				qj = charge[neighbor]
				type_j = typ[neighbor]
				
				dp = pi - pj
				box_func.box_min_dis(dp, self.info.box)
				param = self.params[type_i*self.info.ntypes + type_j]
				fp=[0.0, 0.0]
				rsq=(dp*dp).sum()
#				print(self.func)
				nb_func(rsq, qi, qj, param, fp)
				fij = fp[0]*dp
				result[0] += fij
				result[1][0] += (fij*dp).sum()
				result[1][1] += fp[1]
			forcei = self.info.force[idx]
			virial_potentiali = self.info.virial_potential[idx]
			self.info.force[idx] = forcei+result[0]
			self.info.virial_potential[idx] = virial_potentiali+result[1]
#		print(self.info.virial_potential)

	def setParams(self, type_i, type_j, param, check):	
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)	
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id
		# print(type_i_id, type_j_id, idx1, idx2, self.info.ntypes)
		if isinstance(self.func, str) and self.func == 'lj_coulomb':
			epsilon=param[0]
			sigma=param[1]
			alpha=param[2]
			epsilonr=param[3]
			rcut=param[4]
			p0 = 4.0 * epsilon * math.pow(sigma, int(12))
			p1 = alpha * 4.0 * epsilon * math.pow(sigma, int(6))
			rcutsq = rcut*rcut
			self.params[idx1] = [p0, p1, 138.935/epsilonr, rcutsq]
			self.params[idx2] = [p0, p1, 138.935/epsilonr, rcutsq]
		else:
			self.params[idx1] = param
			self.params[idx2] = param
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True



#--- dpd force
@cuda.jit("void(int32, float32[:, :],float32[:, :], float32[:, :], float32, float32, float32, float32, float32, float32, int32, int32[:], int32[:, :], float32[:, :], float32[:, :], int32, float32, float32, float32)")
def cu_dpd_force(npa, pos, vel, params, box0, box1, box2, box0_half, box1_half, box2_half, ntypes, neighbor_size, neighbor_list, force, virial_potential, seed, T, rsq_dt, one_sixth):
	i = cuda.grid(1)
	if i < npa:
		pix = pos[i][0]
		piy = pos[i][1]
		piz = pos[i][2]	
		
		ti = nb.int32(pos[i][3])
		vi = vel[i]

		result0 = force[i][0]
		result1 = force[i][1]
		result2 = force[i][2]
		result3 = virial_potential[i][0]
		result4 = virial_potential[i][1]	
		for ni in range(nb.int32(0), neighbor_size[i]):
			j = neighbor_list[i][ni]
			dp0 = pix - pos[j][0]
			dp1 = piy - pos[j][1]
			dp2 = piz - pos[j][2]
			tj = nb.int32(pos[j][3])
			# box_func.cu_box_min_dis(dp, box)
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

			pms = params[ti + tj*ntypes]
			alpha = pms[0]
			sigma = pms[1]
			rcutsq = pms[2]
			rcut_inv = pms[3]
			if rsq < rcutsq:
				vj = vel[j]
				dvx = vi[0] - vj[0];
				dvy = vi[1] - vj[1];
				dvz = vi[2] - vj[2];
				
				dot = dp0*dvx +dp1*dvy +dp2*dvz
				rinv = nb.float32(1.0)/math.sqrt(rsq)
				omega = rinv - rcut_inv
				pair_eng = nb.float32(0.25)*alpha*omega*omega*rsq/rcut_inv

				seed1 = nb.int32(0)
				seed2 = nb.int32(0)
				if i < j:
					seed1=i
					seed2=j
				else:
					seed1=j
					seed2=i
				
				rd = potentials.rng.saruprng(nb.int32(seed), nb.int32(seed1), nb.int32(seed2), nb.int32(1), nb.float32(-1.0), nb.float32(1.0))
				gamma = nb.float32(0.5)*sigma*sigma
				Fcfac = alpha*omega
				Fdfac = -dot*gamma*omega*omega
				Frfac = nb.float32(1.7320508)*T*sigma*omega*rd		  
				Forcefac = Fcfac + Fdfac +Frfac*rsq_dt		

				fijx = Forcefac*dp0
				fijy = Forcefac*dp1
				fijz = Forcefac*dp2
				virial = one_sixth * rsq * Forcefac

				result0 += fijx
				result1 += fijy
				result2 += fijz
				result3 += virial
				result4 += pair_eng
		force[i][0] = result0 
		force[i][1] = result1 
		force[i][2] = result2 
		virial_potential[i][0] = result3 
		virial_potential[i][1] = result4

class dpd:
	#define init function
	def __init__(self, info, nlist):
		self.nlist=nlist
		self.info=info
		self.params=[[] for i in range(self.info.ntypes*self.info.ntypes)]
		self.block_size = 256
		self.seed = 12345
		self.T = 1.0
		if self.info.d_velo is None:
			self.info.d_velo = cuda.device_array([self.info.pitch, 3], dtype=np.float32)
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
		
		if self.info.dt<1.0e-7:
			self.rsq_dt = 0.0
		else:
			self.rsq_dt = 1/math.sqrt(self.info.dt)
		
		
		cu_dpd_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.info.d_velo, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
													self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5),self.info.ntypes, d_nlist_size, 
													d_nlist, self.info.d_force, self.info.d_virial_potential, self.seed+timestep, self.T, self.rsq_dt, np.float32(1.0/6.0))
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
	

	def setParams(self, type_i, type_j, alpha, sigma, rcut, check):	
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)	
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id
		self.params[idx1] = [alpha, sigma, rcut*rcut, 1.0/rcut]
		self.params[idx2] = [alpha, sigma, rcut*rcut, 1.0/rcut]

		check[idx1] = True
		check[idx2] = True
		self.params_changed = True

# --- slj force


@cuda.jit("void(int32, float32[:, :],float32[:], float32[:, :], float32, float32, float32, float32, float32, float32, int32, int32[:], int32[:, :], float32[:, :], float32[:, :], float32)")
def cu_slj_force(npa, pos, diameter, params, box0, box1, box2, box0_half, box1_half, box2_half, ntypes, neighbor_size, neighbor_list, force, virial_potential, one_sixth):
	i = cuda.grid(1)
	if i < npa:
		pix = pos[i][0]
		piy = pos[i][1]
		piz = pos[i][2]
		ti = nb.int32(pos[i][3])
		di = diameter[i]

		result0 = force[i][0]
		result1 = force[i][1]
		result2 = force[i][2]
		result3 = virial_potential[i][0]
		result4 = virial_potential[i][1]

		for ni in range(nb.int32(0), neighbor_size[i]):
			j = neighbor_list[i][ni]
			dp0 = pix - pos[j][0]
			dp1 = piy - pos[j][1]
			dp2 = piz - pos[j][2]
			tj = nb.int32(pos[j][3])
			dj = diameter[j]
			
			if dp0 >= box0_half:
				dp0 -= box0
			elif dp0 < -box0_half:
				dp0 += box0

			if dp1 >= box1_half:
				dp1 -= box1
			elif dp1 < -box1_half:
				dp1 += box1

			if dp2 >= box2_half:
				dp2 -= box2
			elif dp2 < -box2_half:
				dp2 += box2

			rsq = dp0*dp0 + dp1*dp1 + dp2*dp2
			pms = params[ti + tj*ntypes]
			# 计算相互作用类型作为索引

			lj1 = pms[0]
			lj2 = pms[1]
			sigma = pms[2]
			rcut = pms[3]
			delta = (di+dj)*nb.float32(0.5)-sigma
			r = math.sqrt(rsq)
			rs = r-delta
			if rs < rcut:
				rinv = nb.float32(1.0)/r
				rsinv = nb.float32(1.0)/rs
				rs2inv = rsinv*rsinv
				rs6inv = rs2inv*rs2inv*rs2inv
				
				force_slj = rinv*rsinv*rs6inv*(nb.float32(12.0)*lj1*rs6inv - nb.float32(6.0)*lj2)
				energy_slj = rs6inv*(lj1*rs6inv-lj2)
				#print(force_slj,energy_slj)
				fijx = force_slj*dp0
				fijy = force_slj*dp1
				fijz = force_slj*dp2

				# virial公式，计算p
				virial = one_sixth * rsq * force_slj

				result0 += fijx
				result1 += fijy
				result2 += fijz
				result3 += virial
				result4 += energy_slj
		force[i][0] = result0
		force[i][1] = result1
		force[i][2] = result2
		virial_potential[i][0] = result3
		virial_potential[i][1] = result4


class slj:
	# define init function
	def __init__(self, info, nlist):
		self.nlist = nlist
		self.info = info
		self.params = [[] for i in range(self.info.ntypes*self.info.ntypes)]
		self.block_size = 256
		self.params_changed = True

	# calculate slj force
	def calculate(self, timestep):
		if self.params_changed:
#			print(self.params)
#			self.info.diameter = self.info.d_diameter.copy_to_host()
#			print(self.info.diameter)
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False

		self.nlist.calculate(timestep)
		d_nlist = self.nlist.data.d_neighbor_list
		d_nlist_size = self.nlist.data.d_neighbor_size
		cu_slj_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.info.d_diameter, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2],
													self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), 
													self.info.box[2]*np.float32(0.5), self.info.ntypes, d_nlist_size,
													d_nlist, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/6.0))
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
	def setParams(self, type_i, type_j, params, check):
		type_i_id = self.info.convert_name_to_id(type_i)
		type_j_id = self.info.convert_name_to_id(type_j)
		idx1 = type_i_id * self.info.ntypes + type_j_id
		idx2 = type_j_id * self.info.ntypes + type_i_id
		self.params[idx1] = params # [epsilon, alpha, sigma, rcut]
		self.params[idx2] = params
		check[idx1] = True
		check[idx2] = True
		self.params_changed = True

