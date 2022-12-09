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

import math
import pygamd.snapshots.box
import numpy as np
import numba as nb
from numba import cuda,njit,jit
from numba.experimental import jitclass


@cuda.jit("void(int32, float32[:, :], int32[:], float32[:], float32[:], int32[:], float32[:, :, :], int32[:])")
def cu_cell_build(npa, pos, dim, box_low_boundary, inv_width, cell_size, cell_list, situation):
	i = cuda.grid(1)
	if i < npa:
		# pi = pos[i]
		pix = pos[i][0]
		piy = pos[i][1]		
		piz = pos[i][2]	
		if math.isnan(pix) or math.isnan(piy) or math.isnan(piz):
			situation[0]=i+nb.int32(1)
			return
			
		if pix<box_low_boundary[0] or pix>=-box_low_boundary[0] or piy<box_low_boundary[1] or piy>=-box_low_boundary[1] or piz<box_low_boundary[2] or piz>=-box_low_boundary[2]:
			situation[1]=i+nb.int32(1)
			return
			
		dpix = pix - box_low_boundary[0]
		dpiy = piy - box_low_boundary[1]
		dpiz = piz - box_low_boundary[2]
		
		ix = nb.int32(dpix*inv_width[0])
		iy = nb.int32(dpiy*inv_width[1])
		iz = nb.int32(dpiz*inv_width[2])		
		if ix == dim[0]:
			ix = nb.int32(0)
		if iy == dim[1]:
			iy = nb.int32(0)			
		if iz == dim[2]:
			iz = nb.int32(0)

		cell_id = iz + dim[2] * (iy + ix * dim[1])
		if cell_id >= cell_list.shape[0]:
			situation[1]=i+nb.int32(1)
			return
		size = cuda.atomic.add(cell_size, cell_id, nb.int32(1))
		if size < cell_list.shape[1]:
			cell_list[cell_id][size][0] = pix
			cell_list[cell_id][size][1] = piy
			cell_list[cell_id][size][2] = piz
			cell_list[cell_id][size][3] = nb.float32(i)
		else:
			cuda.atomic.max(situation, nb.int32(2), size+nb.int32(1))
			
@cuda.jit("void(int32, int32[:])")
def cu_zero(n, array):
	i = cuda.grid(1)
	if i < n:
		array[i] = nb.int32(0)				

@jitclass
class clist:
	#定义构造方法
	def __init__(self, info, rcut):
		self.rcut = rcut
#		print(self.rcut)
		self.info=info
		self.info_box = np.asarray([info.box[0], info.box[1], info.box[2]], dtype=np.float32)
		self.box = np.asarray([info.box[0], info.box[1], info.box[2]], dtype=np.float32)
		self.dim = np.asarray([int(self.box[0]/rcut), int(self.box[1]/rcut), int(self.box[2]/rcut)], dtype=np.int32)
		self.width = np.asarray([self.box[0]/float(self.dim[0]), self.box[1]/float(self.dim[1]), self.box[2]/float(self.dim[2])], dtype=np.float32)
		self.inv_width = np.asarray([1.0/self.width[0], 1.0/self.width[1], 1.0/self.width[2]], dtype=np.float32)
		self.box_low_boundary = np.asarray([-self.box[0]/2.0, -self.box[1]/2.0, -self.box[2]/2.0], dtype=np.float32)
		self.ncell = (self.dim[0]*self.dim[1]*self.dim[2])
		self.cell_adj = np.zeros([self.ncell, 27], dtype = np.int32)
		self.nmax = self.info.npa//self.ncell
		self.nmax = self.nmax + 8 - (self.nmax & 7)
		self.cell_list = np.zeros([self.ncell, self.nmax, 4], dtype = np.float32)
		self.cell_size = np.zeros(self.ncell, dtype = np.int32)
		self.situation = np.zeros(3, dtype = np.int32)
		self.block_size = 64
		
		#--- build cell map
		self.build_cell_map()
		
		#--- device arrays	
		self.d_situation = cuda.to_device(self.situation)	
		self.d_dim = cuda.to_device(self.dim)
		self.d_inv_width = cuda.to_device(self.inv_width)
		self.d_box_low_boundary = cuda.to_device(self.box_low_boundary)
		self.d_cell_size = cuda.to_device(self.cell_size)		
		self.d_cell_list = cuda.to_device(self.cell_list)
		self.d_cell_adj = cuda.to_device(self.cell_adj)	

	def build_cell_map(self):
		for z in range(0, self.dim[2]):
			for y in range(0, self.dim[1]):
				for x in range(0, self.dim[0]):
					cell_id = z + self.dim[2] * (y + x * self.dim[1])
					offset = 0;
					for nz in range(z-1, z+2):
						for ny in range(y-1, y+2):
							for nx in range(x-1, x+2):
								wx = nx % self.dim[0]
								if wx < 0:
									wx += self.dim[0]
									
								wy = ny % self.dim[1]
								if wy < 0:
									wy += self.dim[1]
									
								wz = nz % self.dim[2]
								if wz < 0:
									wz += self.dim[2]									
								
								neigh_cell_id = wz + self.dim[2] * (wy + wx * self.dim[1])
								self.cell_adj[cell_id][offset] = neigh_cell_id
								offset += 1
		for cell_id in range(0, self.ncell):
			self.cell_adj[cell_id]=np.sort(self.cell_adj[cell_id])		
		
	def update_cell_data(self):		
		self.box[0] = self.info.box[0]
		self.box[1] = self.info.box[1]		
		self.box[2] = self.info.box[2]

		dim0 = int(self.box[0]/self.rcut)
		dim1 = int(self.box[1]/self.rcut)			
		dim2 = int(self.box[2]/self.rcut)
		
		if dim0 != self.dim[0] or dim1 != self.dim[1] or dim2 != self.dim[2]:
			self.dim[0] = dim0
			self.dim[1] = dim1
			self.dim[2] = dim2
			
			self.ncell = self.dim[0]*self.dim[1]*self.dim[2]
			self.cell_adj = np.zeros([self.ncell, 27], dtype = np.int32)
			self.cell_list = np.zeros([self.ncell, self.nmax, 4], dtype = np.float32)		
			self.cell_size = np.zeros(self.ncell, dtype = np.int32)
			
			#--- build cell map
			self.build_cell_map()
			#--- device arrays	
			self.d_cell_size = cuda.to_device(self.cell_size)
			self.d_cell_list = cuda.to_device(self.cell_list)
			self.d_cell_adj = cuda.to_device(self.cell_adj)
			self.d_dim = cuda.to_device(self.dim)

		self.width[0] = self.box[0]/float(self.dim[0])
		self.width[1] = self.box[1]/float(self.dim[1])
		self.width[2] = self.box[2]/float(self.dim[2])
	
		self.inv_width[0] = 1.0/self.width[0]
		self.inv_width[1] = 1.0/self.width[1]
		self.inv_width[2] = 1.0/self.width[2]
		
		self.box_low_boundary[0] = -self.box[0]/2.0
		self.box_low_boundary[1] = -self.box[1]/2.0
		self.box_low_boundary[2] = -self.box[2]/2.0	
		#--- device arrays				
		self.d_inv_width = cuda.to_device(self.inv_width)
		self.d_box_low_boundary = cuda.to_device(self.box_low_boundary)
	# @jit # Not possible -> Throw errors
	def calculate(self):
		if (self.info.box[0] != self.box[0] or self.info.box[1] != self.box[1] or self.info.box[2] != self.box[2]) :
			self.update_cell_data()

		build_list = True
		while build_list:
			cu_zero[1, 32](3, self.d_situation)
			cu_zero[math.ceil(self.ncell / self.block_size), self.block_size](self.ncell, self.d_cell_size)

			nblocks = math.ceil(self.info.npa / self.block_size)
			cu_cell_build[nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_dim, self.d_box_low_boundary, self.d_inv_width, self.d_cell_size, self.d_cell_list, self.d_situation)
			self.situation = self.d_situation.copy_to_host()
			if self.situation[0] > 0:
				raise RuntimeError('Error position NaN of particle ')# +str(self.situation[0]-1)+' !')
			if self.situation[1] > 0:
				raise RuntimeError('Error particle ') # +str(self.situation[1]-1)+' out of box !')
			if self.situation[2]>self.nmax:
				self.nmax = self.situation[2]
				self.nmax = self.nmax + 8 - (self.nmax & 7)
				self.cell_list = np.zeros([self.ncell, self.nmax, 4], dtype = np.float32)
				self.d_cell_list = cuda.to_device(self.cell_list)
			else:
				build_list = False
		# self.cell_list = self.d_cell_list.copy_to_host()
		# self.cell_size = self.d_cell_size.copy_to_host()
		# for i in range(0, self.cell_size.shape[0]):
			# if self.cell_size[i]>0:
				# print(self.cell_size[i], self.cell_list[i])
	def speak(self):
		print(self.list)
