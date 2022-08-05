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
import numpy as np
import numba as nb
from numba import cuda
import time

def sort1(output, input):
	output[0] = input[0]
	output[1] = input[3]
	output[2] = input[4]
	output[3] = input[7]
	output[4] = input[6]
	output[5] = input[5]
	output[6] = input[2]
	output[7] = input[1]

def sort2(output, input):
	output[0] = input[0]
	output[1] = input[7]
	output[2] = input[6]
	output[3] = input[1]
	output[4] = input[2]
	output[5] = input[5]
	output[6] = input[4]
	output[7] = input[3]

def sort4(output, input):
	output[0] = input[2]
	output[1] = input[3]
	output[2] = input[0]
	output[3] = input[1]
	output[4] = input[6]
	output[5] = input[7]
	output[6] = input[4]
	output[7] = input[5]

def sort6(output, input):
	output[0] = input[4]
	output[1] = input[3]
	output[2] = input[2]
	output[3] = input[5]
	output[4] = input[6]
	output[5] = input[1]
	output[6] = input[0]
	output[7] = input[7]

def sort8(output, input):
	output[0] = input[6]
	output[1] = input[5]
	output[2] = input[2]
	output[3] = input[1]
	output[4] = input[0]
	output[5] = input[3]
	output[6] = input[4]
	output[7] = input[7]

sort={0: sort1, 1: sort2, 2: sort2, 3: sort4, 4: sort4, 5: sort6, 6: sort6, 7: sort8}
step = np.asarray([[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1], [1, 0, 0]])

# sort rtag
@cuda.jit("void(int32, int32[:], int32[:])")
def cu_sort_rtag(npa, rtag, tag):
	i = cuda.grid(1)
	if i < npa:
		rtag[tag[i]] = i
		
# sort two dimensional int array	
@cuda.jit("void(int32, int32[:], int32[:, :], int32[:, :])")
def cu_sort_int2d(npa, particle_order, temp, array):
	i = cuda.grid(1)
	if i < npa:
		order = particle_order[i]
		for j in range(0, array.shape[1]):
			temp[i][j] = array[order][j]
			
@cuda.jit("void(int32, int32[:, :], int32[:, :])")
def cu_copy_int2d(npa, array, temp):
	i = cuda.grid(1)
	if i < npa:
		for j in range(0, array.shape[1]):
			array[i][j] = temp[i][j]
			
# sort two dimensional float array	
@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :])")
def cu_sort_float2d(npa, particle_order, temp, array):
	i = cuda.grid(1)
	if i < npa:
		order = particle_order[i]
		for j in range(0, array.shape[1]):
			temp[i][j] = array[order][j]
			
@cuda.jit("void(int32, float32[:, :], float32[:, :])")
def cu_copy_float2d(npa, array, temp):
	i = cuda.grid(1)
	if i < npa:
		for j in range(0, array.shape[1]):
			array[i][j] = temp[i][j]

# sort one dimensional int array	
@cuda.jit("void(int32, int32[:], int32[:], int32[:])")
def cu_sort_int1d(npa, particle_order, temp, array):
	i = cuda.grid(1)
	if i < npa:
		order = particle_order[i]
		temp[i] = array[order]
			
@cuda.jit("void(int32, int32[:], int32[:])")
def cu_copy_int1d(npa, array, temp):
	i = cuda.grid(1)
	if i < npa:
		array[i] = temp[i]

# sort one dimensional float array	
@cuda.jit("void(int32, int32[:], float32[:], float32[:])")
def cu_sort_float1d(npa, particle_order, temp, array):
	i = cuda.grid(1)
	if i < npa:
		order = particle_order[i]
		temp[i] = array[order]
			
@cuda.jit("void(int32, float32[:], float32[:])")
def cu_copy_float1d(npa, array, temp):
	i = cuda.grid(1)
	if i < npa:
		array[i] = temp[i]	


@cuda.jit("void(int32, int32, float32[:, :], float32[:], int32[:], int32[:])")
def cu_bin_particle(npa, ngrid, pos, box, order, bin_particle):
	i = cuda.grid(1)
	if i < npa:
		pi = pos[i]			
		fix = (pi[0] + box[0]*nb.float32(0.5))/box[0]
		fiy = (pi[1] + box[1]*nb.float32(0.5))/box[1]		
		fiz = (pi[2] + box[2]*nb.float32(0.5))/box[2]
		
		ix = nb.int32(fix * nb.float32(ngrid))
		iy = nb.int32(fiy * nb.float32(ngrid))
		iz = nb.int32(fiz * nb.float32(ngrid))
		
		if ix == ngrid:
			ix = nb.int32(0)
		if iy == ngrid:
			iy = nb.int32(0)
		if iz == ngrid:
			iz = nb.int32(0)			
		
		bin = ix*(ngrid*ngrid) + iy * ngrid + iz
		bin_particle[i] = order[bin]


class sort_order:
	def __init__(self, info):	
		self.info = info
		if self.info.npa < np.int32(64)*np.int32(64)*np.int32(64):
			self.ngrid = np.int32(64)
		else:
			self.ngrid = np.int32(128)
		self.nelement = self.ngrid*self.ngrid*self.ngrid		
		self.order = np.zeros(self.nelement, dtype=np.int32)
		self.reverse_order = np.zeros(self.nelement, dtype=np.int32)
		self.bin_particle = np.zeros(self.info.npa, dtype=np.int32)
		self.particle_order = np.zeros(self.info.npa, dtype=np.int32)        
		self.count = np.int32(0)
		self.block_size = 256
		self.nblocks = math.ceil(self.info.npa / self.block_size)
		
		cell_order = np.empty(8, dtype=np.int32)
		for i in range(0, 8):
			cell_order[i] = i
		init = np.int32([0, 0, 0])
		self.generate_order(init, self.ngrid, self.ngrid, cell_order, self.reverse_order)
		for i in range(0, self.nelement):
			self.order[self.reverse_order[i]] = i
		# for i in range(0, self.nelement):
			# print(self.order[i])
		self.d_temp_int = cuda.device_array(self.info.npa, dtype=np.int32)
		self.d_temp_float = cuda.device_array(self.info.npa, dtype=np.float32)
		self.d_temp_int3 = cuda.device_array([self.info.npa, 3], dtype=np.int32)
		self.d_temp_float4 = cuda.device_array([self.info.npa, 4], dtype=np.float32)
		self.d_order = cuda.to_device(self.order)
		
		self.d_bin_particle = cuda.device_array(self.info.npa, dtype=np.int32)
		self.d_particle_order = cuda.device_array(self.info.npa, dtype=np.int32)
		
		
	def calculate(self, timestep):
		if self.info.dimension == 3:
			self.get_3d_order()
		self.apply_order()
		

	def get_3d_order(self):
		cu_bin_particle[self.nblocks, self.block_size](self.info.npa, self.ngrid, self.info.d_pos, self.info.d_box, self.d_order, self.d_bin_particle)
		self.bin_particle = self.d_bin_particle.copy_to_host()
		self.particle_order = np.argsort(self.bin_particle)
		self.d_particle_order = cuda.to_device(self.particle_order)

	def generate_order(self, ii, wi, mx, cell_order, reverse_order):
		if wi == np.int32(1):
			reverse_order[self.count] = ii[0]*mx*mx + ii[1]*mx + ii[2]
			self.count += np.int32(1)
		else:
			wi = wi / np.int32(2)
			for m in range(0, 8):
				cell_id = cell_order[m]
				jj = ii + wi * step[cell_id]
				sub_cell_order = np.empty(8, dtype=np.int32)
				sort[m](sub_cell_order, cell_order)
				self.generate_order(jj, wi, mx, sub_cell_order, reverse_order)			
			
	def apply_order(self):
		# sort tag array
		cu_sort_int1d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_int, self.info.d_tag)
		cu_copy_int1d[self.nblocks, self.block_size](self.info.npa, self.info.d_tag, self.d_temp_int)
		
		# sort rtag array		
		cu_sort_rtag[self.nblocks, self.block_size](self.info.npa, self.info.d_rtag, self.info.d_tag)

		# sort position array
		cu_sort_float2d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_float4, self.info.d_pos)		
		cu_copy_float2d[self.nblocks, self.block_size](self.info.npa, self.info.d_pos, self.d_temp_float4)	

		# sort velocity array
		cu_sort_float2d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_float4, self.info.d_vel)		
		cu_copy_float2d[self.nblocks, self.block_size](self.info.npa, self.info.d_vel, self.d_temp_float4)	

		# sort image array
		cu_sort_int2d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_int3, self.info.d_image)		
		cu_copy_int2d[self.nblocks, self.block_size](self.info.npa, self.info.d_image, self.d_temp_int3)
		
		# sort velocity0 array
		if self.info.d_velo is not None:
			cu_sort_float2d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_float4, self.info.d_velo)
			cu_copy_float2d[self.nblocks, self.block_size](self.info.npa, self.info.d_velo, self.d_temp_float4)
			
		# sort diameter array
		if self.info.d_diameter is not None:
			cu_sort_float1d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_float, self.info.d_diameter)
			cu_copy_float1d[self.nblocks, self.block_size](self.info.npa, self.info.d_diameter, self.d_temp_float)

		# sort charge array
		if self.info.d_charge is not None:
			cu_sort_float1d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_float, self.info.d_charge)
			cu_copy_float1d[self.nblocks, self.block_size](self.info.npa, self.info.d_charge, self.d_temp_float)

		# sort body array
		if self.info.d_body is not None:
			cu_sort_int1d[self.nblocks, self.block_size](self.info.npa, self.d_particle_order, self.d_temp_int, self.info.d_body)
			cu_copy_int1d[self.nblocks, self.block_size](self.info.npa, self.info.d_body, self.d_temp_int)				

		for p in self.info.plist:
			p.data.sort()
			p.data.calculate(True)

		for p in self.info.particle_set:
			p.sort()

		if self.info.bond is not None:
			self.info.bond.sort()

		if self.info.angle is not None:
			self.info.angle.sort()

		if self.info.dihedral is not None:
			self.info.dihedral.sort()

			

