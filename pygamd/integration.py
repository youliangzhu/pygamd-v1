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

from pygamd import integrations
from pygamd import chare

class nvt:
	def __init__(self, info, group, method, tau, temperature):
		if method == "nh":
			ci = info.find_comp_info(group)
			if ci is None:
				ci = chare.comp_info(info, group)
				info.comp_info.append(ci)
			self.data=integrations.nvt.nose_hoover(info, ci, tau, temperature)
		self.name="integration"
		self.subname="nvt"

class nve:
	def __init__(self, info, group):
		ci = info.find_comp_info(group)
		if ci is None:
			ci = chare.comp_info(info, group)
			info.comp_info.append(ci)
		self.data = integrations.nve.nve(info, ci)
		self.name = "integration"
		self.subname = "nve"
		
class gwvv:
	def __init__(self, info, group):
		self.data=integrations.dpd.gwvv(info, group)
		self.name="integration"
		self.subname="gwvv"

class bd:
	def __init__(self, info, group, temp):
		self.data=integrations.bd.bd(info, group, temp)
		self.name="integration"
		self.subname="bd"
        
	def setParams(self, typ, gamma):
		self.data.setParams(typ, gamma)	
