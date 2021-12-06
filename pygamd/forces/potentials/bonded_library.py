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

def cu_bond(force_name):
	@cuda.jit(device=True)
	def _harmonic(rsq, param, fp):
		k = param[0]
		r0 = param[1]
		r = math.sqrt(rsq)
		f = k * (r0/r - nb.float32(1.0))
		p = nb.float32(0.5) * k * (r0 - r) * (r0 - r)
		fp[0]=f
		fp[1]=p

	if force_name=="harmonic":
		return _harmonic
		
		
def cu_angle(force_name):
	@cuda.jit(device=True)
	def _harmonic(cos_abc, sin_abc, param, fp):
		k = param[0]
		t0 = param[1]
		dth = math.acos(cos_abc) - t0
		f = k * dth
		p = nb.float32(0.5) * f * dth
		fp[0]=f
		fp[1]=p
		
	@cuda.jit(device=True)
	def _harmonic_cos(cos_abc, sin_abc, param, fp):
		k = param[0]
		t0 = param[1]
		dth = cos_abc - t0
		f = k * dth
		p = nb.float32(0.5) * f * dth
		fp[0]=f
		fp[1]=p		

	if force_name=="harmonic":
		return _harmonic		
	elif force_name=="harmonic_cos":
		return _harmonic_cos
		
def cu_proper(force_name):
	@cuda.jit(device=True)
	def _harmonic(cos_abcd, sin_abcd, param, fp):
		k = param[0]
		cos_phi0 = param[2]
		sin_phi0 = param[3]
		cos_factor = param[4]
		f = cos_factor * (-sin_abcd*cos_phi0 + cos_abcd*sin_phi0)
		p = nb.float32(1.0) + cos_factor * (cos_abcd*cos_phi0 + sin_abcd*sin_phi0)
		fp[0]=-k*f
		fp[1]=k*p
		
	if force_name=="harmonic":
		return _harmonic
		
def cu_improper(force_name):
	@cuda.jit(device=True)
	def _harmonic(cos_abcd, sin_abcd, param, fp):
		k = param[0]
		t0 = param[1];
		dth = math.acos(cos_abcd) - t0
		f = k*dth
		p = nb.float32(0.5) * f * dth
		fp[0]=f
		fp[1]=p

	if force_name=="harmonic":
		return _harmonic

		
