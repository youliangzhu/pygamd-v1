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

from pygamd import tinker
from numba import cuda
import numpy as np
import numba as nb
import time
import math

@cuda.jit("void(int32, float32[:, :], float32[:, :])")
def zero_force_vp(npa, force, virial_potential):
	i = cuda.grid(1)
	if i < npa:
		force[i][0] = nb.float32(0.0)
		force[i][1] = nb.float32(0.0)
		force[i][2] = nb.float32(0.0)
		virial_potential[i][0] = nb.float32(0.0)
		virial_potential[i][1] = nb.float32(0.0)

class dynamics:
	forces=[]
	integrations=[]
	tinkers=[]
	dumps=[]
	def __init__(self, info, dt, sort=True):
		self.info=info
		self.info.dt = dt
		self.sort = sort
		self.start_time = time.perf_counter()
		self.init_timestep = self.info.timestep
		self.last = self.init_timestep
		self.target = self.init_timestep + 200
		self.begin = self.init_timestep		 
		self.end = self.init_timestep
		self.average_tps = 0.0
		self.count_tps = 0
		self.block_size = 64
		if dt>0.0:
			self.sort_period = np.int32(20.0/dt)
		else:
			self.sort_period = np.int32(1000)
		if self.sort:
			self.psort = tinker.psort(info, self.sort_period)

	def run(self, ntimesteps):
		self.end += ntimesteps	
		print("info : --- start to run")
		print("info : from ",self.begin," to ",self.end," timestep")
		
		self.register(self.begin)
		if self.sort:
			self.psort.compute(self.begin)
		for fi in self.forces:
			fi.compute(self.begin)
			
		for di in self.dumps:
			di.compute(self.begin)				
			
		for ts in range(self.begin+1, self.end+1):
			self.register(ts)
			
			for it in self.integrations:
				it.data.firststep(ts)
			
			nblocks = math.ceil(self.info.npa / self.block_size)
			zero_force_vp[nblocks, self.block_size](self.info.npa, self.info.d_force, self.info.d_virial_potential)

			for fi in self.forces:
				fi.compute(ts)
				
			for ii in self.integrations:
				ii.data.secondstep(ts)
				
			for di in self.dumps:
				di.compute(ts)				
				
			for ti in self.tinkers:
				ti.data.compute(ts)
			if self.sort:
				self.psort.compute(ts)
			self.tps_compute(ts)
		
	def add(self, app):
		if app.name=="force":
			self.forces.append(app)

		if app.name=="integration":
			self.integrations.append(app)
			
		if app.name=="tinker":
			self.tinkers.append(app)
			
		if app.name=="dump":
			self.dumps.append(app)			
			
	def remove(self, app):
		if app.name=="force":
			for i in range(0, len(self.forces)):
				if self.forces[i]==app:
					del self.forces[i]
		if app.name=="integration":
			for i in range(0, len(self.integrations)):
				if self.integrations[i]==app:
					del self.integrations[i]
		if app.name=="tinker":
			for i in range(0, len(self.tinkers)):
				if self.tinkers[i]==app:
					del self.tinkers[i]
		if app.name=="dump":
			for i in range(0, len(self.dumps)):
				if self.dumps[i]==app:
					del self.dumps[i]					
					
	def tps_compute(self, timestep):
		if timestep==self.target:	
			self.end_time = time.perf_counter()
			elapsed = self.end_time - self.start_time
			if elapsed <0.00000001 or elapsed>10000000.0:
				raise RuntimeError('Error! tps abnormal with time ', elapsed)			 
			tps = float(self.target-self.last) / elapsed
			all = self.end - timestep
			rseconds = float(all)/tps
			print("info : tps %5.3f"%tps, "  time step %d"%timestep, "  remaining time %5.3f"%rseconds)
			self.average_tps += tps
			self.count_tps += 1;	
			self.last = self.target
			self.target += int(20.0*tps)+1
			self.start_time = self.end_time
			if self.target>self.end and self.count_tps!=0:
				print("info : average tps %5.3f"%(self.average_tps/float(self.count_tps)))

	def register(self, timestep):
		for key in self.info.compute_properties:
			self.info.compute_properties[key] = False
			
		for ii in self.integrations:
			ii.data.register(timestep)
				
		for di in self.dumps:
			di.register(timestep)				
 
					
