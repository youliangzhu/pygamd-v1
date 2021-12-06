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

from pygamd import forces
from pygamd import plist

class nonbonded:
	def __init__(self, info, rcut, func, exclusion=None):
		nl = info.find_plist(rcut, exclusion)
		if nl is None:
			nl = plist.neighbor(info, rcut, exclusion)
			info.plist.append(nl)	
		self.info = info
		self.data=forces.pair.pair(info, nl, func)
		self.name="force"
		self.subname="nonbonded"
		self.check=[False for i in range(self.info.ntypes*self.info.ntypes)]
		self.first_compute = True
		
	def setParams(self, type_i, type_j, param):
		self.data.setParams(type_i, type_j, param, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)
	
	def check_parameters(self):
		for i in range(self.info.ntypes):
			for j in range(i, self.info.ntypes):
				idx = i * self.info.ntypes + j
				if not self.check[idx]:
					raise RuntimeError('Error! the parameters between type ',self.info.typemap[i],' and type ',self.info.typemap[j],' has not been set!')
class dpd:
	def __init__(self, info, rcut=1.0):
		nl = info.find_plist(rcut, None)
		if nl is None:
			nl = plist.neighbor(info, rcut, None)
			info.plist.append(nl)
		self.info = info
		self.data=forces.pair.dpd(info, nl)
		self.name="force"
		self.subname="dpd"
		self.rcut = rcut
		self.check=[False for i in range(self.info.ntypes*self.info.ntypes)]
		self.first_compute = True		
		
	def setParams(self, type_i, type_j, alpha, sigma):
		self.data.setParams(type_i, type_j, alpha, sigma, self.rcut, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False	
		self.data.calculate(timestep)

	def check_parameters(self):
		for i in range(self.info.ntypes):
			for j in range(i, self.info.ntypes):
				idx = i * self.info.ntypes + j
				if not self.check[idx]:
					raise RuntimeError('Error! the parameters between type ',self.info.typemap[i],' and type ',self.info.typemap[j],' have not been set!')
		
class bond:
	def __init__(self, info, func):	
		self.info = info
		self.data=forces.bond.bond(info, func)
		self.name="force"
		self.subname="bond"
		self.check=[False for i in range(self.info.bond.nbtypes)]
		self.first_compute = True
		
	def setParams(self, bond_type, param):
		self.data.setParams(bond_type, param, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)
	
	def check_parameters(self):
		for i in range(self.info.bond.nbtypes):
			if not self.check[i]:
				raise RuntimeError('Error! the parameters for bond type ',self.info.bond.typemap[i],' have not been set!')		
				
				
class angle:
	def __init__(self, info, func):	
		self.info = info
		self.data=forces.angle.angle(info, func)
		self.name="force"
		self.subname="angle"
		self.check=[False for i in range(self.info.angle.natypes)]
		self.first_compute = True
		
	def setParams(self, angle_type, param):
		self.data.setParams(angle_type, param, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)
	
	def check_parameters(self):
		for i in range(self.info.angle.natypes):
			if not self.check[i]:
				raise RuntimeError('Error! the parameters for angle type ',self.info.angle.typemap[i],' have not been set!')		
				
				
class dihedral:
	def __init__(self, info, func):	
		self.info = info
		self.data=forces.dihedral.dihedral(info, func)
		self.name="force"
		self.subname="dihedral"
		self.check=[False for i in range(self.info.dihedral.ndtypes)]
		self.first_compute = True
		
	def setParams(self, dihedral_type, param, term='proper'):
		self.data.setParams(dihedral_type, param, term, self.check)
		
	def setCosFactor(factor):
		self.data.setCosFactor(factor)
		
	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)
	
	def check_parameters(self):
		for i in range(self.info.dihedral.ndtypes):
			if not self.check[i]:
				raise RuntimeError('Error! the parameters for dihedral type ',self.info.dihedral.typemap[i],' have not been set!')		
				
				
				
				
		