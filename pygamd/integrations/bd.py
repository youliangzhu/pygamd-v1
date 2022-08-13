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
from pygamd import chare
import pygamd.snapshots.box as box_func
import numpy as np
import numba as nb
from numba import cuda
import math 

@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :], float32[:, :], int32[:, :], float32[:], float32, float32)")
def cu_first_step(nme, member, pos, vel, force, image, box, dt, dtsq):
	i = cuda.grid(1)
	if i < nme:
		idx = member[i]
		pi = pos[idx]
		vi = vel[idx]
		ai = force[idx]
		mi = vi[3]
		ii = image[idx]

		ai[0] /= mi
		ai[1] /= mi
		ai[2] /= mi	

		pi[0] += vi[0]*dt + nb.float32(0.5)*ai[0]*dtsq
		pi[1] += vi[1]*dt + nb.float32(0.5)*ai[1]*dtsq
		pi[2] += vi[2]*dt + nb.float32(0.5)*ai[2]*dtsq

		vi[0] += nb.float32(0.5)*ai[0]*dt
		vi[1] += nb.float32(0.5)*ai[1]*dt
		vi[2] += nb.float32(0.5)*ai[2]*dt

		box_func.cu_box_wrap(pi, box, ii)
		pos[idx][0] = pi[0]
		pos[idx][1] = pi[1]
		pos[idx][2] = pi[2]
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]
		
		image[idx][0] = ii[0]
		image[idx][1] = ii[1]
		image[idx][2] = ii[2]

@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :], float32[:, :], float32[:], int32, float32, int32, float32, float32)")
def cu_second_step(nme, member, pos, vel, force, params, seed, T, D, dt, dt_inv):
	i = cuda.grid(1)
	if i < member.shape[0]:
		idx = member[i]
		vi = vel[idx]
		fi = force[idx]
		mi = vi[3]
		ti = nb.int32(pos[idx][3])
		gamma = params[ti]
		coeff = math.sqrt(nb.float32(6.0) * gamma * T * dt_inv)

		rdx = potentials.rng.saruprng(nb.int32(seed), nb.int32(idx), nb.int32(seed+1), nb.int32(1), nb.float32(-1.0), nb.float32(1.0))
		rdy = potentials.rng.saruprng(nb.int32(seed), nb.int32(seed+2), nb.int32(idx), nb.int32(1), nb.float32(-1.0), nb.float32(1.0))
		rdz = potentials.rng.saruprng(nb.int32(seed+3), nb.int32(idx), nb.int32(seed), nb.int32(1), nb.float32(-1.0), nb.float32(1.0))
		
		rfx = rdx*coeff - gamma*vi[0]
		rfy = rdy*coeff - gamma*vi[1]
		rfz = rdz*coeff - gamma*vi[2]
		
		if D == 2:
			rfz = nb.float32(0.0)
			
		fi[0] += rfx 
		fi[1] += rfy
		fi[2] += rfz		
		
		m_inv = nb.float32(1.0) / mi
		
		vi[0] = vi[0] + nb.float32(0.5)*fi[0]*m_inv*dt
		vi[1] = vi[1] + nb.float32(0.5)*fi[1]*m_inv*dt
		vi[2] = vi[2] + nb.float32(0.5)*fi[2]*m_inv*dt
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]
		
		force[idx][0] = fi[0]
		force[idx][1] = fi[1]
		force[idx][2] = fi[2]

class bd:
	#define init function
	def __init__(self, info, group, temp):
		self.info=info
		self.block_size=64
		self.group = group
		self.temp = temp
		self.seed = nb.int32(12345)
		self.params = [1.0 for i in range(self.info.ntypes)]
		self.params_changed = True
		
		self.ps = info.find_particle_set(group)
		if self.ps is None:
			self.ps = chare.particle_set(info, group)
			info.particle_set.append(self.ps)
			
	def setParams(self, typ, gamma):
		typ_id = self.info.convert_name_to_id(typ)
		self.params[typ_id] = gamma
		self.params_changed = True
		
	#calculate non-bonded force
	def firststep(self, timestep):
		nblocks = math.ceil(self.ps.nme / self.block_size)
		cu_first_step[nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_pos, self.info.d_vel, self.info.d_force, 
												self.info.d_image, self.info.d_box, self.info.dt, self.info.dt*self.info.dt)

	def secondstep(self, timestep):
		if self.params_changed:
			self.h_params = np.asarray(self.params, dtype=np.float32)
			self.d_params = cuda.to_device(self.h_params)
			self.params_changed = False
		dt_inv = nb.float32(0.0)
		if self.info.dt > nb.float32(0.0):
			dt_inv = nb.float32(1.0)/self.info.dt
		nblocks = math.ceil(self.ps.nme / self.block_size)
		
		cu_second_step[nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_pos, self.info.d_vel, self.info.d_force, 
												self.d_params, self.seed + timestep, self.temp, self.info.dimension, self.info.dt, dt_inv)
	def register(self, timestep):
		return


