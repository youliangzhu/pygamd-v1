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
import os

@cuda.jit("void(int32, int32[:], float64[:, :], float32[:, :], float64[:, :], int32[:, :], float32[:], float32, float32[:, :])")
def cu_first_step(nme, member, pos, vel, force, image, box, dt, testai):
	i = cuda.grid(1)
	if i < nme:
		idx = member[i]
		pi = pos[idx]
		vi = vel[idx]
		mi = vi[3]
		ai = force[idx]
		ii = image[idx]
		ai[0] = force[idx][0] /mi
		ai[1] = force[idx][1] /mi
		ai[2] = force[idx][2] /mi

		testai[idx][0] = force[idx][0]
		testai[idx][1] = force[idx][1]
		testai[idx][2] = force[idx][2]

		vi[0] = vi[0] + nb.float32(0.5)*ai[0]*dt
		vi[1] = vi[1] + nb.float32(0.5)*ai[1]*dt
		vi[2] = vi[2] + nb.float32(0.5)*ai[2]*dt
		
		pi[0] += vi[0]*dt
		pi[1] += vi[1]*dt
		pi[2] += vi[2]*dt
		
		#box_func.cu_box_wrap(pi, box, ii)
		if pi[0] >=box[0]*nb.float32(0.5):
			pi[0] -= box[0]
			ii[0] += nb.int32(1)
		elif pi[0]<-box[0]*nb.float32(0.5):
			pi[0] += box[0]
			ii[0] -= nb.int32(1)
		
		if pi[1]>=box[1]*nb.float32(0.5):
			pi[1] -= box[1]
			ii[1] += nb.int32(1)		
		elif pi[1]<-box[1]*nb.float32(0.5):
			pi[1] += box[1]
			ii[1] -= nb.int32(1)		

		if pi[2]>=box[2]*nb.float32(0.5):
			pi[2] -= box[2]
			ii[2] += nb.int32(1)		
		elif pi[2]<-box[2]*nb.float32(0.5):
			pi[2] += box[2]
			ii[2] -= nb.int32(1)		

		pos[idx][0] = pi[0]
		pos[idx][1] = pi[1]
		pos[idx][2] = pi[2]
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]
		
		image[idx][0] = ii[0]
		image[idx][1] = ii[1]
		image[idx][2] = ii[2]

@cuda.jit("void(int32, int32[:], float32[:, :], float64[:, :], float32)")
def cu_second_step(nme, member, vel, force, dt):
	i = cuda.grid(1)
	if i < member.shape[0]:
		idx = member[i]
		vi = vel[idx]
		ai = force[idx]
		mi = vi[3]
		
		ai[0] = ai[0] / mi
		ai[1] = ai[1] / mi
		ai[2] = ai[2] / mi
		
		vi[0] = vi[0] + nb.float32(0.5)*ai[0]*dt
		vi[1] = vi[1] + nb.float32(0.5)*ai[1]*dt
		vi[2] = vi[2] + nb.float32(0.5)*ai[2]*dt
		
		vel[idx][0] = vi[0]
		vel[idx][1] = vi[1]
		vel[idx][2] = vi[2]
			
class nve:
	def __init__(self, info, ci):
		self.info=info
		self.ci=ci
		self.block_size=64
		self.xi=0.0
		self.ai = self.info.force
		self.d_ai = cuda.to_device(self.ai)
		self.info.d_pos = cuda.to_device(nb.float64(self.info.pos))
	def firststep(self, timestep):
		nblocks = math.ceil(self.ci.ps.nme / self.block_size)

		cu_first_step[nblocks, self.block_size](self.ci.ps.nme, self.ci.ps.d_member, self.info.d_pos, self.info.d_vel, self.info.d_force, 
												self.info.d_image, self.info.d_box, self.info.dt, self.d_ai)

	def secondstep(self, timestep):
		self.ci.calculate(timestep)
		nblocks = math.ceil(self.ci.ps.nme / self.block_size)

		cu_second_step[nblocks, self.block_size](self.ci.ps.nme, self.ci.ps.d_member, self.info.d_vel, self.info.d_force, self.info.dt)

	def register(self, timestep):
		return
