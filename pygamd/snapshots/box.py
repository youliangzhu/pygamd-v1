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
import numba as nb

# box functions
def box_wrap(pos, box, image):
	if pos[0] >=box [0]/2.0:
		pos[0] -= box[0]
		image[0] += 1
	elif pos[0]<-box[0]/2.0:
		pos[0] += box[0]
		image[0] -= 1
		
	if pos[1]>=box[1]/2.0:
		pos[1] -= box[1]
		image[1] += 1		
	elif pos[1]<-box[1]/2.0:
		pos[1] += box[1]
		image[1] -= 1		

	if pos[2]>=box[2]/2.0:
		pos[2] -= box[2]
		image[2] += 1		
	elif pos[2]<-box[2]/2.0:
		pos[2] += box[2]
		image[2] -= 1		

def box_min_dis(dpos, box):
	if dpos[0]>=box[0]/2.0:
		dpos[0] -= box[0]
	elif dpos[0]<-box[0]/2.0:
		dpos[0] += box[0]
		
	if dpos[1]>=box[1]/2.0:
		dpos[1] -= box[1]
	elif dpos[1]<-box[1]/2.0:
		dpos[1] += box[1]

	if dpos[2]>=box[2]/2.0:
		dpos[2] -= box[2]
	elif dpos[2]<-box[2]/2.0:
		dpos[2] += box[2]

@cuda.jit("void(float32[:], float32[:], int32[:])", device=True)		
def cu_box_wrap(pos, box, image):
	if pos[0] >=box [0]*nb.float32(0.5):
		pos[0] -= box[0]
		image[0] += nb.int32(1)
	elif pos[0]<-box[0]*nb.float32(0.5):
		pos[0] += box[0]
		image[0] -= nb.int32(1)
		
	if pos[1]>=box[1]*nb.float32(0.5):
		pos[1] -= box[1]
		image[1] += nb.int32(1)		
	elif pos[1]<-box[1]*nb.float32(0.5):
		pos[1] += box[1]
		image[1] -= nb.int32(1)		

	if pos[2]>=box[2]*nb.float32(0.5):
		pos[2] -= box[2]
		image[2] += nb.int32(1)		
	elif pos[2]<-box[2]*nb.float32(0.5):
		pos[2] += box[2]
		image[2] -= nb.int32(1)			

@cuda.jit("void(float32[:], float32[:])", device=True)
def cu_box_min_dis(dpos, box):
	if dpos[0]>=box[0]*nb.float32(0.5):
		dpos[0] -= box[0]
	elif dpos[0]<-box[0]*nb.float32(0.5):
		dpos[0] += box[0]
		
	if dpos[1]>=box[1]*nb.float32(0.5):
		dpos[1] -= box[1]
	elif dpos[1]<-box[1]*nb.float32(0.5):
		dpos[1] += box[1]

	if dpos[2]>=box[2]*nb.float32(0.5):
		dpos[2] -= box[2]
	elif dpos[2]<-box[2]*nb.float32(0.5):
		dpos[2] += box[2]		
		
		
		
		