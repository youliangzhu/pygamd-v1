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
from numba import cuda
import numpy as np

@cuda.jit(device=True)
def saruprng(seed1, seed2, seed3, nbound, low, high):
	seed3^=(seed1<<7)^(seed2>>6)
	seed2+=(seed1>>4)^(seed3>>15)
	seed1^=(seed2<<9)+(seed3<<8)
	seed3^=0xA5366B4D*((seed2>>11) ^ (seed1<<1))
	seed2+=0x72BE1579*((seed1<<4)  ^ (seed3>>16))
	seed1^=0X3F38A6ED*((seed3>>5)  ^ (((np.int32)(seed2))>>22))
	seed2+=seed1*seed3
	seed1+=seed3 ^ (seed2>>2)
	seed2^=((np.int32)(seed2))>>17
	
	statex  = 0x79dedea3*(seed1^(((np.int32)(seed1))>>14))
	statey = (statex+seed2) ^ (((np.int32)(statex))>>8)
	statex  = statex + (statey*(statey^0xdddf97f5))
	statey = 0xABCB96F7 + (statey>>1)
	sum = 0
	delta=0x9e3779b9
	k0 = 0xA341316C
	k1 = 0xC8013EA4
	k2 = 0xAD90777D
	k3 = 0x7E95761E
	k4 = 0x6957f5a7
	for i in range(0, nbound):
		sum += delta;
		statex += ((statey<<4)+k0)^(statey+sum)^((statey>>5)+k1)
		statey += ((statex<<4)+k2)^(statex+sum)^((statex>>5)+k3)
	v=np.uint32((statex ^ (statex>>26))+statey)
	u32 = (v^(v>>20))*k4
	TWO_N32 = 0.232830643653869628906250e-9
	return np.uint32(u32)*(TWO_N32*(high-low))+low


		
