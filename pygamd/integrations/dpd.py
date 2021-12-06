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

@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :], float32[:, :], float32[:, :], int32[:, :], float32[:], float32, float32)")
def cu_first_step(nme, member, pos, vel, velo, force, image, box, lam, dt):
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

		vox = vi[0] + lam * dt * ai[0]
		voy = vi[1] + lam * dt * ai[1]
		voz = vi[2] + lam * dt * ai[2]		

		vi[0] = vi[0] + nb.float32(0.5)*ai[0]*dt
		vi[1] = vi[1] + nb.float32(0.5)*ai[1]*dt
		vi[2] = vi[2] + nb.float32(0.5)*ai[2]*dt

		pi[0] += vi[0]*dt
		pi[1] += vi[1]*dt				
		pi[2] += vi[2]*dt

		box_func.cu_box_wrap(pi, box, ii)
		pos[idx][0] = pi[0]
		pos[idx][1] = pi[1]
		pos[idx][2] = pi[2]
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]
		
		velo[idx][0] = vox
		velo[idx][1] = voy
		velo[idx][2] = voz		
		
		image[idx][0] = ii[0]
		image[idx][1] = ii[1]
		image[idx][2] = ii[2]

@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :], float32)")
def cu_second_step(nme, member, vel, force, dt):
	i = cuda.grid(1)
	if i < member.shape[0]:
		idx = member[i]
		vi = vel[idx]
		ai = force[idx]
		mi = vi[3]
		
		ai[0] /= mi
		ai[1] /= mi
		ai[2] /= mi
		
		vi[0] = vi[0] + nb.float32(0.5)*ai[0]*dt
		vi[1] = vi[1] + nb.float32(0.5)*ai[1]*dt
		vi[2] = vi[2] + nb.float32(0.5)*ai[2]*dt
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]


class gwvv:
	#define init function
	def __init__(self, info, group):
		self.info=info
		self.block_size=64
		self.group = group
		
		self.ps = info.find_particle_set(group)
		if self.ps is None:
			self.ps = chare.particle_set(info, group)
			info.particle_set.append(self.ps)		

		self.lam = 0.65
		if self.info.d_velo is None:
			self.info.d_velo = cuda.device_array([self.info.pitch, 3], dtype=np.float32)	
		
	#calculate non-bonded force
	def firststep(self, timestep):
		nblocks = math.ceil(self.ps.nme / self.block_size)
		cu_first_step[nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_pos, self.info.d_vel, self.info.d_velo, 
												self.info.d_force, self.info.d_image, self.info.d_box, self.lam, self.info.dt)
		
	def secondstep(self, timestep):
		nblocks = math.ceil(self.ps.nme / self.block_size)
		cu_second_step[nblocks, self.block_size](self.ps.nme, self.ps.d_member, self.info.d_vel, self.info.d_force, self.info.dt)
		
	def register(self, timestep):
		return


