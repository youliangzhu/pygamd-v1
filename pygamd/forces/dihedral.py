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

try:
	deviceFunction = cuda.compiler.DeviceFunctionTemplate
except:
	deviceFunction = cuda.compiler.DeviceDispatcher

minimum_value = nb.float32(0.001)

def dihedral_force(cu_func_pro, cu_func_imp):
	if isinstance(cu_func_pro, str):
		cu_func_pro = potentials.bonded_library.cu_proper(cu_func_pro)
	if isinstance(cu_func_imp, str):
		cu_func_imp = potentials.bonded_library.cu_improper(cu_func_imp)			
	if not isinstance(cu_func_pro, deviceFunction) or not isinstance(cu_func_imp, deviceFunction):
		raise RuntimeError('Error dihedral_force device function!')
	@cuda.jit("void(int32, float32[:, :], float32[:, :], float32, float32, float32, float32, float32, float32, int32[:], int32[:, :, :], float32[:, :], float32[:, :], float32, float32, int32)")
	def cu_dihedral_force(npa, pos, params, box0, box1, box2, box0_half, box1_half, box2_half, dihedral_size, dihedral_list, force, virial_potential, one_four, one_twelve, param_len):
		i = cuda.grid(1)
		if i < npa:
			pi = pos[i]
			# result = cuda.local.array(5, dtype = nb.float32)
			result0 = force[i][0]
			result1 = force[i][1]
			result2 = force[i][2]
			result3 = virial_potential[i][0]
			result4 = virial_potential[i][1]
			
			for di in range(nb.int32(0), dihedral_size[i]):
				j = dihedral_list[i][di][0]
				k = dihedral_list[i][di][1]
				l = dihedral_list[i][di][2]				
				type  = dihedral_list[i][di][3]
				order = dihedral_list[i][di][4]
				
				pj = pos[j]
				pk = pos[k]
				pl = pos[l]
				
				if order == 0:
					pa = pi
					pb = pj
					pc = pk
					pd = pl	
				elif order == 1:
					pa = pj
					pb = pi
					pc = pk
					pd = pl
				elif order == 2:
					pa = pj
					pb = pk
					pc = pi
					pd = pl
				elif order == 3:
					pa = pj
					pb = pk
					pc = pl
					pd = pi
					

				d_ab0 = pa[0] - pb[0]
				d_ab1 = pa[1] - pb[1]
				d_ab2 = pa[2] - pb[2]
				
				d_cb0 = pc[0] - pb[0]
				d_cb1 = pc[1] - pb[1]
				d_cb2 = pc[2] - pb[2]

        
				d_dc0 = pd[0] - pc[0]
				d_dc1 = pd[1] - pc[1]
				d_dc2 = pd[2] - pc[2]
		

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


				if d_dc0>=box0_half:
					d_dc0 -= box0
				elif d_dc0<-box0_half:
					d_dc0 += box0
					
				if d_dc1>=box1_half:
					d_dc1 -= box1
				elif d_dc1<-box1_half:
					d_dc1 += box1
				
				if d_dc2>=box2_half:
					d_dc2 -= box2
				elif d_dc2<-box2_half:
					d_dc2 += box2

				pms = params[type]
				term = pms[param_len-1]
				fp = cuda.local.array(2, dtype = nb.float32)
				
				if term ==0:    # proper
					d_cbm0 = -d_cb0
					d_cbm1 = -d_cb1
					d_cbm2 = -d_cb2

					aa0 = d_ab1*d_cbm2 - d_ab2*d_cbm1
					aa1 = d_ab2*d_cbm0 - d_ab0*d_cbm2
					aa2 = d_ab0*d_cbm1 - d_ab1*d_cbm0
					
					bb0 = d_dc1*d_cbm2 - d_dc2*d_cbm1;
					bb1 = d_dc2*d_cbm0 - d_dc0*d_cbm2;
					bb2 = d_dc0*d_cbm1 - d_dc1*d_cbm0;
					
					aa_rsq = aa0*aa0 + aa1*aa1 + aa2*aa2
					bb_rsq = bb0*bb0 + bb1*bb1 + bb2*bb2
					g_rsq = d_cbm0*d_cbm0 + d_cbm1*d_cbm1 + d_cbm2*d_cbm2
					g_r = math.sqrt(g_rsq)
					
					g_rinv = nb.float32(0.0)
					aa_rsqinv = nb.float32(0.0)
					bb_rsqinv = nb.float32(0.0)
					
					if g_r > nb.float32(0.0):
						g_rinv = nb.float32(1.0)/g_r
						
					if aa_rsq > nb.float32(0.0): 
						aa_rsqinv = nb.float32(1.0)/aa_rsq
						
					if bb_rsq > nb.float32(0.0):
						bb_rsqinv = nb.float32(1.0)/bb_rsq
						
					ab_rinv = math.sqrt(aa_rsqinv*bb_rsqinv)
					cos_abcd = (aa0*bb0 + aa1*bb1 + aa2*bb2)*ab_rinv
					sin_abcd = g_r*ab_rinv*(aa0*d_dc0 + aa1*d_dc1 + aa2*d_dc2)
					# print(sin_abcd, g_r, ab_rinv, aa0*d_dc0 + aa1*d_dc1 + aa2*d_dc2)
					
					if cos_abcd > nb.float32(1.0):
						cos_abcd = nb.float32(1.0)
						
					if cos_abcd < -nb.float32(1.0): 
						cos_abcd = -nb.float32(1.0)

					cu_func_pro(cos_abcd, sin_abcd, pms, fp)
					f = fp[0]
					
					fg = d_ab0*d_cbm0 + d_ab1*d_cbm1 + d_ab2*d_cbm2
					hg = d_dc0*d_cbm0 + d_dc1*d_cbm1 + d_dc2*d_cbm2
					
					fga = fg*aa_rsqinv*g_rinv
					hgb = hg*bb_rsqinv*g_rinv
					gaa = -aa_rsqinv*g_r
					gbb = bb_rsqinv*g_r
					
					df0 = gaa*aa0
					df1 = gaa*aa1
					df2 = gaa*aa2
					
					dg0 = fga*aa0 - hgb*bb0
					dg1 = fga*aa1 - hgb*bb1
					dg2 = fga*aa2 - hgb*bb2
					
					dth0 = gbb*bb0
					dth1 = gbb*bb1
					dth2 = gbb*bb2

					ss0 = f*dg0
					ss1 = f*dg1
					ss2 = f*dg2
					
					ffa0 = f*df0
					ffa1 = f*df1
					ffa2 = f*df2
					
					ffb0 = ss0 - ffa0
					ffb1 = ss1 - ffa1
					ffb2 = ss2 - ffa2
					
					ffd0 = f*dth0
					ffd1 = f*dth1
					ffd2 = f*dth2
					
					ffc0 = -ss0 - ffd0
					ffc1 = -ss1 - ffd1
					ffc2 = -ss2 - ffd2
					
				elif term ==1:   # improper
					r1 = np.float32(1.0)/math.sqrt(d_ab0*d_ab0 + d_ab1*d_ab1 + d_ab2*d_ab2)
					r2 = np.float32(1.0)/math.sqrt(d_cb0*d_cb0 + d_cb1*d_cb1 + d_cb2*d_cb2)
					r3 = np.float32(1.0)/math.sqrt(d_dc0*d_dc0 + d_dc1*d_dc1 + d_dc2*d_dc2)
					
					ss1 = r1 * r1
					ss2 = r2 * r2
					ss3 = r3 * r3
					
					# cosine and sin of the angle between the planes
					c0 =  (d_ab0*d_dc0 + d_ab1*d_dc1 + d_ab2*d_dc2) * r1 * r3
					c1 =  (d_ab0*d_cb0 + d_ab1*d_cb1 + d_ab2*d_cb2) * r1 * r2
					c2 = -(d_dc0*d_cb0 + d_dc1*d_cb1 + d_dc2*d_cb2) * r3 * r2
					
					s1 = np.float32(1.0) - c1*c1
					if s1 < minimum_value:
						s1 = minimum_value
					s1 = np.float32(1.0) / s1
					
					s2 = np.float32(1.0) - c2*c2
					if s2 < minimum_value:
						s2 = minimum_value
					s2 = np.float32(1.0) / s2
					
					s12 = math.sqrt(s1*s2)
					cos_abcd = (c1*c2 + c0) * s12;
					
					if cos_abcd > np.float32(1.0):
						cos_abcd = np.float32(1.0)
						
					if cos_abcd < -np.float32(1.0):
						cos_abcd = -np.float32(1.0)
					
					sin_abcd = math.sqrt(np.float32(1.0) - cos_abcd*cos_abcd)
					if sin_abcd < minimum_value:
						sin_abcd = minimum_value


					cu_func_imp(cos_abcd, sin_abcd, pms, fp)
					f = fp[0]
					
					f = -f /sin_abcd
					cos_abcd = cos_abcd * f
					s12 = s12 * f
					a11 = cos_abcd*ss1*s1;
					a22 = -ss2 * (np.float32(2.0)*c0*s12 - cos_abcd*(s1+s2))
					a33 = cos_abcd*ss3*s2;
					
					a12 = -r1*r2*(c1*cos_abcd*s1 + c2*s12)
					a13 = -r1*r3*s12;
					a23 = r2*r3*(c2*cos_abcd*s2 + c1*s12)
					
					ss0 = a22*d_cb0 + a23*d_dc0 + a12*d_ab0
					ss1 = a22*d_cb1 + a23*d_dc1 + a12*d_ab1
					ss2 = a22*d_cb2 + a23*d_dc2 + a12*d_ab2
					
					# calculate the forces for each particle
					ffa0 = a12*d_cb0 + a13*d_dc0 + a11*d_ab0
					ffa1 = a12*d_cb1 + a13*d_dc1 + a11*d_ab1
					ffa2 = a12*d_cb2 + a13*d_dc2 + a11*d_ab2
					
					ffb0 = -ss0 - ffa0
					ffb1 = -ss1 - ffa1
					ffb2 = -ss2 - ffa2
					
					ffd0 = a23*d_cb0 + a33*d_dc0 + a13*d_ab0
					ffd1 = a23*d_cb1 + a33*d_dc1 + a13*d_ab1
					ffd2 = a23*d_cb2 + a33*d_dc2 + a13*d_ab2
					
					ffc0 = ss0 - ffd0
					ffc1 = ss1 - ffd1
					ffc2 = ss2 - ffd2	
				
				if order == 0:
					result0 += ffa0
					result1 += ffa1
					result2 += ffa2
					
				if order == 1:
					result0 += ffb0
					result1 += ffb1
					result2 += ffb2
					
				if order == 2:
					result0 += ffc0
					result1 += ffc1
					result2 += ffc2
					
				if order == 3:
					result0 += ffd0
					result1 += ffd1
					result2 += ffd2

				vx = d_ab0*ffa0 + d_cb0*ffc0 + (d_dc0 + d_cb0)*ffd0
				vy = d_ab1*ffa1 + d_cb1*ffc1 + (d_dc1 + d_cb1)*ffd1
				vz = d_ab2*ffa2 + d_cb2*ffc2 + (d_dc2 + d_cb2)*ffd2
				virial = one_twelve*(vx + vy + vz)
				
				potential = fp[1] * one_four
				result3 += virial
				result4 += potential
			force[i][0] = result0 
			force[i][1] = result1 
			force[i][2] = result2 
			virial_potential[i][0] = result3 
			virial_potential[i][1] = result4
	return cu_dihedral_force

class dihedral:
	#define init function
	def __init__(self, info, func):
		self.func=func
		self.info=info
		self.params=[[] for i in range(self.info.dihedral.ndtypes)]
		self.block_size = 256
		self.dihedral_force = dihedral_force(func, func)
		self.cos_factor = -1.0
		self.params_changed = True
	#calculate dihedral force
	def calculate(self, timestep):
		if self.params_changed:
			# update cos_factor for harmonic
			if isinstance(self.func, str) and self.func == 'harmonic':
				npar = len(self.params[0])
				for i in range(self.info.dihedral.ndtypes):
					self.params[i][npar-2] = self.cos_factor
		
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.nblocks = math.ceil(self.info.npa / self.block_size)
			self.param_len = self.h_params.shape[1]
			self.params_changed = False
		# for i in range(0, self.info.npa):
			# if self.info.dihedral.dihedral_size[i] >0:
				# print(i, self.info.dihedral.dihedral_size[i], self.info.dihedral.dihedral_list[i])
		# start = time.time()	
		self.dihedral_force[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_params, self.info.box[0], self.info.box[1], self.info.box[2], 
														self.info.box[0]*np.float32(0.5), self.info.box[1]*np.float32(0.5), self.info.box[2]*np.float32(0.5), self.info.dihedral.d_dihedral_size,  
														self.info.dihedral.d_dihedral_list, self.info.d_force, self.info.d_virial_potential, np.float32(1.0/4.0), np.float32(1.0/12.0), self.param_len)
		# cuda.synchronize()
		# end1 = time.time()
		# print("force",end1 - start)											   
		# self.info.force = self.info.d_force.copy_to_host()
		# self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
		# for i in range(0, self.info.force.shape[0]):
			# print(i, self.info.force[i], self.info.virial_potential[i])
			
	def setCosFactor(factor):
		self.cos_factor = factor
		self.params_changed = True

	def setParams(self, dihedral_type, param, term, check):
		type_id = self.info.dihedral.convert_name_to_id(dihedral_type)
		if term == "proper":
			term_id = 0
		elif term =="improper":
			term_id = 1
		else:
			raise RuntimeError("Error, dihedral term ", term, " the candidates are 'dihedral' and 'improper' !")
		if isinstance(self.func, str) and self.func == 'harmonic':
			if len(param) != 2:
				raise RuntimeError("Error, the number of harmonic parameters is not equal to 2!")
			k=param[0]
			t0=param[1]
			t0_rad = math.pi*t0/180.0
			self.params[type_id] = [k, t0_rad, math.cos(t0_rad), math.sin(t0_rad), self.cos_factor, term_id]		
		elif isinstance(self.func, deviceFunction):
			param.append(term_id)
			self.params[type_id] = param

		check[type_id] = True
		self.params_changed = True
