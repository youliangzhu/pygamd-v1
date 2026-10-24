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

from bisect import bisect_right


class Variant:
	def __init__(self):
		self._offset = 0

	def getValue(self, timestep):
		return 0.0

	def setOffset(self, offset):
		if offset < 0:
			raise ValueError("Variant offset must be non-negative")
		self._offset = int(offset)


class VariantConst(Variant):
	def __init__(self, val):
		super().__init__()
		self._val = float(val)

	def getValue(self, timestep):
		return self._val


class VariantLinear(Variant):
	def __init__(self):
		super().__init__()
		self._times = []
		self._values = []
		self._last_index = 0

	def setPoint(self, timestep, val):
		timestep = int(timestep)
		val = float(val)

		if timestep < 0:
			raise ValueError("VariantLinear timestep must be non-negative")

		idx = bisect_right(self._times, timestep)

		# Overwrite existing point at the same timestep.
		if idx > 0 and self._times[idx - 1] == timestep:
			self._values[idx - 1] = val
			self._last_index = 0
			return

		self._times.insert(idx, timestep)
		self._values.insert(idx, val)
		self._last_index = 0

	def getValue(self, timestep):
		timestep = int(timestep)

		# Match GALAMOST Variant offset semantics.
		if timestep < self._offset:
			timestep = 0
		else:
			timestep -= self._offset

		if len(self._times) == 0:
			raise RuntimeError("Error getting variant value: no points specified to VariantLinear")

		if len(self._times) == 1:
			return self._values[0]

		if timestep < self._times[0]:
			return self._values[0]

		if timestep >= self._times[-1]:
			return self._values[-1]

		# Cache the interpolation bracket for monotonically increasing timesteps.
		i = self._last_index
		if not (0 <= i < len(self._times) - 1 and self._times[i] <= timestep < self._times[i + 1]):
			i = bisect_right(self._times, timestep) - 1
			self._last_index = i

		ta = self._times[i]
		tb = self._times[i + 1]
		va = self._values[i]
		vb = self._values[i + 1]

		f = float(timestep - ta) / float(tb - ta)
		return (1.0 - f) * va + f * vb

	def __call__(self, timestep):
		return self.getValue(timestep)
