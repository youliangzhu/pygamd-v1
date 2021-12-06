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

import numpy as np
import numba as nb
from numba import cuda
import math 

#temperature
@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :],  float32[:], int32)")
def cu_sums1(nme, member, vel, virial_potential, coll, nblocks):
	sm = cuda.shared.array(256, nb.float32)
	i = cuda.grid(1)
	tx = cuda.threadIdx.x
	temp = nb.float32(0.0)
	if i < nme:
		idx = member[i]
		vi = vel[idx]
		mi = vi[3]
		temp = mi*(vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2])
		
	sm[tx] = temp
	cuda.syncthreads()
	
	offs = cuda.blockDim.x >> nb.int32(1);
	while offs > nb.int32(0):
		if tx < offs:
			sm[tx] += sm[tx + offs]		
		offs >>= nb.int32(1)
		cuda.syncthreads()
		
	if tx == nb.int32(0):
		coll[cuda.blockIdx.x] = sm[0]

@cuda.jit("void(float32[:], float32[:], int32)")	
def cu_final_sums1(coll, result, nblocks):
	sm = cuda.shared.array(512, nb.float32)
	tx = cuda.threadIdx.x	
	final_sum_temp = nb.float32(0.0)
	for offset in range(nb.int32(0), nblocks, cuda.blockDim.x):
		cuda.syncthreads()
		if offset + tx < nblocks:
			temp = coll[offset + tx]
			sm[tx] = temp	
		else:
			sm[tx] = nb.float32(0.0)
		cuda.syncthreads()

		offs = cuda.blockDim.x >> nb.int32(1)
		while offs > nb.int32(0):
			if tx < offs:
				sm[tx] += sm[tx + offs]				
			offs >>= nb.int32(1)
			cuda.syncthreads()

		if tx == nb.int32(0):
			final_sum_temp += sm[0]

	if tx == nb.int32(0):
		result[0] = final_sum_temp


# temperature and pressure
@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :],  float32[:], int32)")
def cu_sums2(nme, member, vel, virial_potential, coll, nblocks):
	block_size=256
	sm = cuda.shared.array(256*2, nb.float32)
	i = cuda.grid(1)
	tx = cuda.threadIdx.x
	temp = nb.float32(0.0)
	virial = nb.float32(0.0)
	if i < nme:
		idx = member[i]
		vi = vel[idx]
		mi = vi[3]
		vp = virial_potential[idx]
		temp = mi*(vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2])
		virial = vp[0]
		
	sm[tx] = temp
	sm[tx+block_size] = virial
	cuda.syncthreads()
	
	offs = cuda.blockDim.x >> 1;
	while offs > 0:
		if tx < offs:
			sm[tx] += sm[tx + offs]
			sm[tx + block_size] += sm[tx + block_size + offs]			
		offs >>= 1
		cuda.syncthreads()
		
	if tx == 0:
		coll[cuda.blockIdx.x] = sm[0]
		coll[cuda.blockIdx.x + nblocks] = sm[block_size]

@cuda.jit("void(float32[:], float32[:], int32)")	
def cu_final_sums2(coll, result, nblocks):
	block_size=512
	sm = cuda.shared.array(512*2, nb.float32)
	tx = cuda.threadIdx.x	
	final_sum_temp = nb.float32(0.0)
	final_sum_virial = nb.float32(0.0)
	for offset in range(0, nblocks, cuda.blockDim.x):
		cuda.syncthreads()
		if offset + tx < nblocks:
			temp = coll[offset + tx]
			virial = coll[offset + nblocks+ tx]

			sm[tx] = temp
			sm[tx + block_size] = virial		
		else:
			sm[tx] = nb.float32(0.0)
			sm[tx + block_size] = nb.float32(0.0)
		cuda.syncthreads()

		offs = cuda.blockDim.x >> 1
		while offs > 0:
			if tx < offs:
				sm[tx] += sm[tx + offs]
				sm[tx + block_size] += sm[tx + block_size + offs]				
			offs >>= 1
			cuda.syncthreads()

		if tx == 0:
			final_sum_temp += sm[0]
			final_sum_virial += sm[block_size]

	if tx == 0:
		result[0] = final_sum_temp
		result[1] = final_sum_virial


# momentum, temperature, pressure and potential
@cuda.jit("void(int32, int32[:], float32[:, :], float32[:, :],  float32[:], int32)")
def cu_sums4(nme, member, vel, virial_potential, coll, nblocks):
	block_size=256
	block_size2=block_size*2
	block_size3=block_size*3
	block_size4=block_size*4
	block_size5=block_size*5
	
	sm = cuda.shared.array(256*6, nb.float32)
	i = cuda.grid(1)
	tx = cuda.threadIdx.x
	
	temp = nb.float32(0.0)
	virial = nb.float32(0.0)
	potential = nb.float32(0.0)
	mx = nb.float32(0.0)
	my = nb.float32(0.0)
	mz = nb.float32(0.0)
	
	if i < nme:
		idx = member[i]
		vi = vel[idx]
		mi = vi[3]
		vp = virial_potential[idx]
		temp = mi*(vi[0]*vi[0] + vi[1]*vi[1] + vi[2]*vi[2])
		virial = vp[0]
		potential = vp[1]
		mx = vi[0]*mi
		my = vi[1]*mi
		mz = vi[2]*mi		
		
	sm[tx] = temp
	sm[tx+block_size] = virial
	sm[tx+block_size2] = potential	
	sm[tx+block_size3] = mx	
	sm[tx+block_size4] = my	
	sm[tx+block_size5] = mz		
	cuda.syncthreads()
	
	offs = cuda.blockDim.x >> 1;
	while offs > 0:
		if tx < offs:
			sm[tx] += sm[tx + offs]
			sm[tx + block_size] += sm[tx + block_size + offs]
			sm[tx + block_size2] += sm[tx + block_size2 + offs]
			sm[tx + block_size3] += sm[tx + block_size3 + offs]
			sm[tx + block_size4] += sm[tx + block_size4 + offs]
			sm[tx + block_size5] += sm[tx + block_size5 + offs]			
		offs >>= 1
		cuda.syncthreads()
		
	if tx == 0:
		coll[cuda.blockIdx.x] = sm[0]
		coll[cuda.blockIdx.x + nblocks] = sm[block_size]
		coll[cuda.blockIdx.x + nblocks*2] = sm[block_size2]
		coll[cuda.blockIdx.x + nblocks*3] = sm[block_size3]
		coll[cuda.blockIdx.x + nblocks*4] = sm[block_size4]
		coll[cuda.blockIdx.x + nblocks*5] = sm[block_size5]
	
@cuda.jit("void(float32[:], float32[:], int32)")	
def cu_final_sums4(coll, result, nblocks):
	block_size=512
	block_size2=block_size*2
	block_size3=block_size*3
	block_size4=block_size*4
	block_size5=block_size*5	
	sm = cuda.shared.array(512*6, nb.float32)
	tx = cuda.threadIdx.x
	
	final_sum_temp = nb.float32(0.0)
	final_sum_virial = nb.float32(0.0)
	final_sum_potential = nb.float32(0.0)
	final_sum_mx = nb.float32(0.0)
	final_sum_my = nb.float32(0.0)
	final_sum_mz = nb.float32(0.0)	
	for offset in range(0, nblocks, cuda.blockDim.x):
		cuda.syncthreads()
		if offset + tx < nblocks:
			temp = coll[offset + tx]
			virial = coll[offset + nblocks+ tx]
			potential = coll[offset + 2*nblocks + tx]
			mx = coll[offset + 3*nblocks + tx]
			my = coll[offset + 4*nblocks + tx]
			mz = coll[offset + 5*nblocks + tx]
			
			sm[tx] = temp
			sm[tx + block_size] = virial
			sm[tx + block_size2] = potential
			sm[tx + block_size3] = mx
			sm[tx + block_size4] = my
			sm[tx + block_size5] = mz			
		else:
			sm[tx] = nb.float32(0.0)
			sm[tx + block_size] = nb.float32(0.0)
			sm[tx + block_size2] = nb.float32(0.0)
			sm[tx + block_size3] = nb.float32(0.0)
			sm[tx + block_size4] = nb.float32(0.0)
			sm[tx + block_size5] = nb.float32(0.0)			
		cuda.syncthreads()

		offs = cuda.blockDim.x >> 1
		while offs > 0:
			if tx < offs:
				sm[tx] += sm[tx + offs]
				sm[tx + block_size] += sm[tx + block_size + offs]
				sm[tx + block_size2] += sm[tx + block_size2 + offs]
				sm[tx + block_size3] += sm[tx + block_size3 + offs]
				sm[tx + block_size4] += sm[tx + block_size4 + offs]
				sm[tx + block_size5] += sm[tx + block_size5 + offs]				
			offs >>= 1
			cuda.syncthreads()

		if tx == 0:
			final_sum_temp += sm[0]
			final_sum_virial += sm[block_size]
			final_sum_potential += sm[block_size2]
			final_sum_mx += sm[block_size3]
			final_sum_my += sm[block_size4]
			final_sum_mz += sm[block_size5]			
	if tx == 0:
		result[0] = final_sum_temp
		result[1] = final_sum_virial
		result[2] = final_sum_potential			
		result[3] = final_sum_mx
		result[4] = final_sum_my
		result[5] = final_sum_mz
		
