'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
'''

from pygamd import tinkers


class psort:
	def __init__(self, info, period):
		self.data = tinkers.hilbert.sort_order(info)
		self.name = "tinker"
		self.subname = "sort"
		self.period = int(period)

	def setGrid(self, mx, my, mz):
		self.data.setGrid(mx, my, mz)

	def compute(self, timestep):
		if self.period > 0 and timestep % self.period == 0:
			self.data.calculate(timestep)
