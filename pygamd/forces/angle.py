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

minimum_value = nb.float32(0.001)

def angle_force(cu_func):
	if isinstance(cu_func, str):
		cu_func = potentials.bonded_library.cu_angle(cu_func)
	if not isinstance(cu_func, cuda.compiler.DeviceDispatcher):
		raise RuntimeError('Error angle_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:, :], float32, float32, float32, float32, float32, float32, int32[:], int32[:, :, :], float32[:, :], float32[:, :], float32, float32)")
	def cu_angle_force(npa, pos, params, box0, box1, box2, box0_half, box1_half, box2_half, angle_size, angle_list, force, virial_potential, one_three, one_sixth):
		i = cuda.grid(1)
		if i < npa:
			pi = pos[i]
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]	
			for ai in range(nb.int32(0), angle_size[i]):
				j = angle_list[i][ai][0]
				k = angle_list[i][ai][1]
				type  = angle_list[i][ai][2]
				order = angle_list[i][ai][3]
				
				pj = pos[j]
				pk = pos[k]
				
				if order == 0:
					pa = pi
					pb = pj
					pc = pk
					
				if order == 1:
					pa = pj
					pb = pi
					pc = pk
					
				if order == 2:
					pa = pj
					pb = pk
					pc = pi					

				d_ab0 = pa[0] - pb[0]
				d_ab1 = pa[1] - pb[1]
				d_ab2 = pa[2] - pb[2]
				
				d_cb0 = pc[0] - pb[0]
				d_cb1 = pc[1] - pb[1]
				d_cb2 = pc[2] - pb[2]

				if d_ab0>=box0_half:
					d_ab0 -= box0
				elif d_ab0<-box0_half:
					d_ab0 += box0
					
				if d_ab1>=box1_half:
					d_ab1 -= box1
				elif d_ab1<-box1_half:
					d_ab1 += box1
				
				if d_ab2>=box2_half:
					d_ab2 -= box2
				elif d_ab2<-box2_half:
					d_ab2 += box2	

				
				if d_cb0>=box0_half:
					d_cb0 -= box0
				elif d_cb0<-box0_half:
					d_cb0 += box0
					
				if d_cb1>=box1_half:
					d_cb1 -= box1
				elif d_cb1<-box1_half:
					d_cb1 += box1
				
				if d_cb2>=box2_half:
					d_cb2 -= box2
				elif d_cb2<-box2_half:
					d_cb2 += box2						
					
				rsq_ab = d_ab0*d_ab0 + d_ab1*d_ab1 + d_ab2*d_ab2
				r_ab = math.sqrt(rsq_ab)
				rsq_cb = d_cb0*d_cb0 + d_cb1*d_cb1 + d_cb2*d_cb2				
				r_cb = math.sqrt(rsq_cb)
				
				cos_abc = d_ab0*d_cb0 + d_ab1*d_cb1 + d_ab2*d_cb2
				cos_abc /= r_ab*r_cb
				
				if cos_abc > nb.float32(1.0):
					cos_abc = nb.float32(1.0)
				if cos_abc < -nb.float32(1.0): 
					cos_abc = -nb.float32(1.0)
					
					
				sin_abc = math.sqrt(nb.float32(1.0) - cos_abc*cos_abc)
				if sin_abc < minimum_value: 
					sin_abc = minimum_value
				sin_abc = nb.float32(1.0)/sin_abc					
				
				pms = params[type]
				fp = cuda.local.array(2, dtype = nb.float32)
				cu_func(cos_abc, sin_abc, pms, fp)

				a = -fp[0] * sin_abc
				a11 = a*cos_abc/rsq_ab
				a12 = -a / (r_ab*r_cb)
				a22 = a*cos_abc / rsq_cb
        
				fab0 = a11*d_ab0 + a12*d_cb0
				fab1 = a11*d_ab1 + a12*d_cb1
				fab2 = a11*d_ab2 + a12*d_cb2
				
				fcb0 = a22*d_cb0 + a12*d_ab0
				fcb1 = a22*d_cb1 + a12*d_ab1
				fcb2 = a22*d_cb2 + a12*d_ab2
				
				if order == 0:
					result0 += fab0
					result1 += fab1
					result2 += fab2
					
				if order == 1:
					result0 -= fab0 + fcb0
					result1 -= fab1 + fcb1
					result2 -= fab2 + fcb2
					
				if order == 2:
					result0 += fcb0
					result1 += fcb1
					result2 += fcb2
				
				vx = d_ab0*fab0 + d_cb0*fcb0
				vy = d_ab1*fab1 + d_cb1*fcb1
				vz = d_ab2*fab2 + d_cb2*fcb2	
				virial = one_sixth*(vx + vy + vz)
				# if i==35 and ai == 0:
					# print(i, ai, vy, a, cos_abc, rsq_cb)
				potential = fp[1] * one_three
				result3 += virial
				result4 += potential
			force[i][0] = result0 
			force[i][1] = result1 
			force[i][2] = result2 
			virial_potential[i][0] = result3 
			virial_potential[i][1] = result4
	return cu_angle_force

class angle:
	#define init function
	def __init__(self, info, func):
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.angle.natypes)]
		self.block_size = 256
		self.angle_force = angle_force(func)
		self.params_changed = True
	#calculate angle force
	def calculate(self, timestep):
		if self.params_changed:
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.params_changed = False

		# for i in range(0, self.info.npa):
			# if self.info.angle.angle_size[i] >0:
				# print(i, self.info.angle.angle_size[i], self.info.angle.angle_list[i])
		# start = time.time()	
		self.angle_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
														self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.angle.d_angle_size,  
														self.info.angle.d_angle_list, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/3.0), np.float32(1.0/6.0))
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)											   
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])

	def setParams(self, angle_type, param, check):
		type_id = self.info.angle.convert_name_to_id(angle_type)
		if isinstance(self.func, str) and self.func == 'harmonic':
			if len(param) != 2:
				raise RuntimeError('Error, the number of harmonic parameters is not equal to 2!')
			k=param[0]
			t0=param[1]
			t0_rad = math.pi*t0/180.0
			self.params[type_id] = [k, t0_rad]
		elif isinstance(self.func, str) and self.func == 'harmonic_cos':
			if len(param) != 2:
				raise RuntimeError('Error, the number of harmonic parameters is not equal to 2!')
			k=param[0]
			t0=param[1]
			t0_rad = math.pi*t0/180.0
			self.params[type_id] = [k, math.cos(t0_rad)]			
		elif isinstance(self.func, cuda.compiler.DeviceFunctionTemplate):
			self.params[type_id] = param

		check[type_id] = True
		self.params_changed = True
