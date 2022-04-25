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

# from pygamd import dumps
from pygamd import chare
from numba import cuda
import numpy as np
import math

class data:
	def __init__(self, info, group, file, period):
		self.name="dump"
		self.subname="data"
		self.period = period
		self.info = info
		self.ci = info.find_comp_info(group)
		if self.ci is None:
			self.ci = chare.comp_info(info, group)
			info.comp_info.append(self.ci)
		self.file = file
		with open(self.file, 'w+') as f:
			f.write("timestep\tmomentum\ttemperature\tpressure\tpotential\n")	
		
	def compute(self, timestep):
		if timestep%self.period==0:
	
			self.ci.calculate(timestep)
			# output data
			with open(self.file, 'a+') as f:
				#f.write(str(timestep)+"  "+str(self.ci.momentum)+"	 "+str(self.ci.temp)+"	"+str(self.ci.pressure)+"  "+str(self.ci.potential)+"\n")
				f.write('%10d\t%.6e\t%.6f\t%.6e\t%.6f\n' % (timestep, self.ci.momentum, self.ci.temp, self.ci.pressure, self.ci.potential))
				
	def register(self, timestep):
		if timestep%self.period==0:
			self.info.compute_properties['momentum']=True
			self.info.compute_properties['temperature']=True			
			self.info.compute_properties['pressure']=True			
			self.info.compute_properties['potential']=True			
			

class mst:
	def __init__(self, info, group, file, period, properties=None, split=False):
		self.name="dump"
		self.subname="mst"
		self.period = period
		self.info = info
		self.file = file

		if self.file is None:
			raise RuntimeError('Error! file name is not given!')
			
		file_list = self.file.split('.')
		if file_list[len(file_list)-1] != 'mst':
			self.file += '.mst'

		self.group = group
		self.ps = info.find_particle_set(group)
		if self.ps is None:
			self.ps = chare.particle_set(info, group)
			info.particle_set.append(self.ps)
		
		self.properties=properties
		self.split = split
		self.frame = 0
		self.first_compute=True
		self.dump = {'position':True, 'type':True, 'mass':True, 'image':True, 'box':True,
					'velocity':False, 'force':False, 'potential':False, 'virial':False,
					'charge':False, 'body':False, 'diameter':False, 'init':False, 'cris':False,
					'bond':True, 'angle':True, 'dihedral':True}
		
		if self.properties is not None:
			for key in self.dump:
				self.dump[key]=False
				
			for p in self.properties:
				self.dump[p] = True

	def compute(self, timestep):
		if timestep%self.period==0:
			# output data
			if not self.split:
				if self.first_compute:
					with open(self.file, 'w') as f:
						f.write("mst_version 1.0 #Copyright You-Liang Zhu\n")
						f.write("invariant_data\n")
						self.write_properties(f, 0, False)
						f.write("variant_data\n")
					self.first_compute = False

				with open(self.file, 'a') as f:
					f.write('frame\t%d\n' % (self.frame))
					self.write_properties(f, timestep, True)
					f.write('frame_end\n')
			else:
				ip = self.file.find('.mst')
				if ip == -1:
					raise RuntimeError('Error! file name!', self.file)
				file_new = self.file[0:ip]+'.'+str(timestep).rjust(10,'0')+'.mst'
				with open(file_new, 'w') as f:
					f.write("mst_version 1.0 #Copyright You-Liang Zhu\n")
					self.write_properties(f, timestep, False)
					self.write_properties(f, timestep, True)
					f.write('mst_end\n')
			self.frame += 1
			
	def write_properties(self, f, timestep, variant_indicator):
		if variant_indicator:
			f.write('\ttimestep\n')	
			f.write('\t\t%d\n' % (timestep))
		else:
			f.write('\tnum_particles\n')
			f.write('\t\t%d\n' % (self.ps.nme))
			f.write('\tdimension\n')
			f.write('\t\t%d\n' % (self.info.dimension))
			
		self.info.rtag = self.info.d_rtag.copy_to_host()
		if self.info.variant["box"]==variant_indicator and self.dump["box"]:
			f.write('\tbox\n')
			f.write('\t\t%10.6f\t%10.6f\t%10.6f\n' % (self.info.box[0], self.info.box[1], self.info.box[2]))
			
		if self.info.variant["position"]==variant_indicator and self.dump["position"]:
			self.info.pos = self.info.d_pos.copy_to_host()
			f.write('\tposition\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\t%10.6f\t%10.6f\n' % (self.info.pos[idx][0], self.info.pos[idx][1], self.info.pos[idx][2]))
				
		if self.info.variant["type"]==variant_indicator and self.dump["type"]:
			self.info.pos = self.info.d_pos.copy_to_host()
			f.write('\ttype\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				typid = np.int32(self.info.pos[idx][3])
				f.write('\t\t'+self.info.typemap[typid]+'\n')
				
		if self.info.variant["velocity"]==variant_indicator and self.dump["velocity"]:
			self.info.vel = self.info.d_vel.copy_to_host()
			f.write('\tvelocity\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\t%10.6f\t%10.6f\n' % (self.info.vel[idx][0], self.info.vel[idx][1], self.info.vel[idx][2]))
				
		if self.info.variant["image"]==variant_indicator and self.dump["image"]:
			self.info.image = self.info.d_image.copy_to_host()
			f.write('\timage\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10d\t%10d\t%10d\n' % (self.info.image[idx][0], self.info.image[idx][1], self.info.image[idx][2]))
				
		if self.info.variant["mass"]==variant_indicator and self.dump["mass"]:
			self.info.vel = self.info.d_vel.copy_to_host()
			f.write('\tmass\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\n' % (self.info.vel[idx][3]))
				
		if self.info.variant["force"]==variant_indicator and self.dump["force"]:
			self.info.force = self.info.d_force.copy_to_host()
			f.write('\tforce\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\t%10.6f\t%10.6f\n' % (self.info.force[idx][0], self.info.force[idx][1], self.info.force[idx][2]))
				
		if self.info.variant["virial"]==variant_indicator and self.dump["virial"]:
			self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
			f.write('\tvirial\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\n' % (self.info.virial_potential[idx][0]))
				
		if self.info.variant["potential"]==variant_indicator and self.dump["potential"]:
			self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
			f.write('\tpotential\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\n' % (self.info.virial_potential[idx][1]))
				
		if self.info.charge is not None and self.info.variant["charge"]==variant_indicator and self.dump["charge"]:
			self.info.charge = self.info.d_charge.copy_to_host()
			f.write('\tcharge\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\n' % (self.info.charge[idx]))
				
		if self.info.body is not None and self.info.variant["body"]==variant_indicator and self.dump["body"]:
			self.info.body = self.info.d_body.copy_to_host()
			f.write('\tbody\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10d\n' % (self.info.body[idx]))

		if self.info.diameter is not None and self.info.variant["diameter"]==variant_indicator and self.dump["diameter"]:
			self.info.diameter = self.info.d_diameter.copy_to_host()
			f.write('\tdiameter\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10.6f\n' % (self.info.diameter[idx]))

		if self.info.init is not None and self.info.variant["init"]==variant_indicator and self.dump["init"]:
			self.info.init = self.info.d_init.copy_to_host()
			f.write('\tinit\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10d\n' % (self.info.init[idx]))

		if self.info.cris is not None and self.info.variant["cris"]==variant_indicator and self.dump["cris"]:
			self.info.cris = self.info.d_cris.copy_to_host()
			f.write('\tcris\n')
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('\t\t%10d\n' % (self.info.cris[idx]))
				
		if self.info.bond is not None and self.info.variant["bond"]==variant_indicator and self.dump["bond"]:
			f.write('\tbond\n')
			for i in range(0, len(self.info.bond.bonds)):
				bi = self.info.bond.bonds[i]
				f.write('\t\t%s\t%10d\t%10d\n' % (bi[0], bi[1], bi[2]))
				
		if self.info.angle is not None and self.info.variant["angle"]==variant_indicator and self.dump["angle"]:
			f.write('\tangle\n')
			for i in range(0, len(self.info.angle.angles)):
				ai = self.info.angle.angles[i]
				f.write('\t\t%s\t%10d\t%10d\t%10d\n' % (ai[0], ai[1], ai[2], ai[3]))
				
		if self.info.dihedral is not None and self.info.variant["dihedral"]==variant_indicator and self.dump["dihedral"]:
			f.write('\tdihedral\n')
			for i in range(0, len(self.info.dihedral.dihedrals)):
				di = self.info.dihedral.dihedrals[i]
				f.write('\t\t%s\t%10d\t%10d\t%10d\t%10d\n' % (di[0], di[1], di[2], di[3], di[4]))
				
	def register(self, timestep):
		return

class xml:
	def __init__(self, info, group, file, period, properties=None, split=True):
		self.name="dump"
		self.subname="xml"
		self.period = period
		self.info = info
		self.file = file

		if self.file is None:
			raise RuntimeError('Error! file name is not given!')
			
		file_list = self.file.split('.')
		if file_list[len(file_list)-1] != 'xml':
			self.file += '.xml'

		self.group = group
		self.ps = info.find_particle_set(group)
		if self.ps is None:
			self.ps = chare.particle_set(info, group)
			info.particle_set.append(self.ps)		
		
		self.properties=properties
		self.split = split
		self.frame = 0
		self.first_compute=True
		self.dump = {'position':True, 'type':True, 'mass':True, 'image':True,
					'velocity':False, 'force':False, 'potential':False, 'virial':False,
					'charge':False, 'body':False, 'diameter':False, 'init':False, 'cris':False,
					'bond':True, 'angle':True, 'dihedral':True}
		
		if self.properties is not None:
			for key in self.dump:
				self.dump[key]=False

			for p in self.properties:
				self.dump[p] = True

	def compute(self, timestep):
		if timestep%self.period==0:
			# output data
			if not self.split:
				if self.first_compute:
					with open(self.file, 'w') as f:
						f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
						f.write('<galamost_xml version="1.3">\n')
					self.first_compute = False
				with open(self.file, 'a') as f:
					self.write_properties(f, timestep)
			else:
				ip = self.file.find('.xml')
				if ip == -1:
					raise RuntimeError('Error! file name!', self.file)
				file_new = self.file[0:ip]+'.'+str(timestep).rjust(10,'0')+'.xml'
				with open(file_new, 'w') as f:
					f.write('<?xml version="1.0" encoding="UTF-8"?>\n')
					f.write('<galamost_xml version="1.3">\n')
					self.write_properties(f, timestep)
					f.write('</galamost_xml>\n')
			self.frame += 1
			
	def write_properties(self, f, timestep):
		f.write('<configuration time_step="%d" dimensions="%d" natoms="%d" >\n' % (timestep, self.info.dimension, self.ps.nme))
		f.write('<box lx="%10.6f" ly="%10.6f" lz="%10.6f"/>\n' % (self.info.box[0], self.info.box[1], self.info.box[2]))
		self.info.rtag = self.info.d_rtag.copy_to_host()
		if self.dump["position"]:
			self.info.pos = self.info.d_pos.copy_to_host()
			f.write('<position num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\t%10.6f\t%10.6f\n' % (self.info.pos[idx][0], self.info.pos[idx][1], self.info.pos[idx][2]))
			f.write('</position>\n')
			
		if self.dump["type"]:
			self.info.pos = self.info.d_pos.copy_to_host()
			f.write('<type num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				typid = np.int32(self.info.pos[idx][3])
				f.write(self.info.typemap[typid]+'\n')
			f.write('</type>\n')
			
		if self.dump["velocity"]:
			self.info.vel = self.info.d_vel.copy_to_host()
			f.write('<velocity num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\t%10.6f\t%10.6f\n' % (self.info.vel[idx][0], self.info.vel[idx][1], self.info.vel[idx][2]))
			f.write('</velocity>\n')
			
		if self.dump["image"]:
			self.info.image = self.info.d_image.copy_to_host()
			f.write('<image num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10d\t%10d\t%10d\n' % (self.info.image[idx][0], self.info.image[idx][1], self.info.image[idx][2]))	
			f.write('</image>\n')
			
		if self.dump["mass"]:
			self.info.vel = self.info.d_vel.copy_to_host()
			f.write('<mass num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\n' % (self.info.vel[idx][3]))
			f.write('</mass>\n') 
			
		if self.dump["force"]:
			self.info.force = self.info.d_force.copy_to_host()
			f.write('<force num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\t%10.6f\t%10.6f\n' % (self.info.force[idx][0], self.info.force[idx][1], self.info.force[idx][2]))
			f.write('</force>\n') 
			
		if self.dump["virial"]:
			self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
			f.write('<virial num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\n' % (self.info.virial_potential[idx][0]))
			f.write('</virial>\n') 
			
		if self.dump["potential"]:
			self.info.virial_potential = self.info.d_virial_potential.copy_to_host()
			f.write('<potential num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\n' % (self.info.virial_potential[idx][1]))
			f.write('</potential>\n')
			
		if self.info.charge is not None and self.dump["charge"]:
			self.info.charge = self.info.d_charge.copy_to_host()
			f.write('<charge num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\n' % (self.info.charge[idx]))
			f.write('</charge>\n')
			
		if self.info.body is not None and self.dump["body"]:
			self.info.body = self.info.d_body.copy_to_host()
			f.write('<body num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10d\n' % (self.info.body[idx]))
			f.write('</body>\n')

		if self.info.diameter is not None and self.dump["diameter"]:
			self.info.diameter = self.info.d_diameter.copy_to_host()
			f.write('<diameter num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10.6f\n' % (self.info.diameter[idx]))
			f.write('</diameter>\n')

		if self.info.init is not None and self.dump["init"]:
			self.info.init = self.info.d_init.copy_to_host()
			f.write('<init num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10d\n' % (self.info.init[idx]))
			f.write('</init>\n')

		if self.info.cris is not None and self.dump["cris"]:
			self.info.cris = self.info.d_cris.copy_to_host()
			f.write('<cris num="%d">\n' % (self.ps.nme))
			for i in range(0, self.ps.nme):
				tag = self.ps.member[i]
				idx = self.info.rtag[tag] 
				f.write('%10d\n' % (self.info.cris[idx]))
			f.write('</cris>\n')			
			
		if self.info.bond is not None and self.dump["bond"]:
			f.write('<bond num="%d">\n' % (len(self.info.bond.bonds)))
			for i in range(0, len(self.info.bond.bonds)):
				bi = self.info.bond.bonds[i]
				f.write('%s\t%10d\t%10d\n' % (bi[0], bi[1], bi[2]))
			f.write('</bond>\n')
			
		if self.info.angle is not None and self.dump["angle"]:
			f.write('<angle num="%d">\n' % (len(self.info.angle.angles)))
			for i in range(0, len(self.info.angle.angles)):
				ai = self.info.angle.angles[i]
				f.write('%s\t%10d\t%10d\t%10d\n' % (ai[0], ai[1], ai[2], ai[3]))
			f.write('</angle>\n')
			
		if self.info.dihedral is not None and self.dump["dihedral"]:
			f.write('<dihedral num="%d">\n' % (len(self.info.dihedral.dihedrals)))
			for i in range(0, len(self.info.dihedral.dihedrals)):
				di = self.info.dihedral.dihedrals[i]
				f.write('%s\t%10d\t%10d\t%10d\t%10d\n' % (di[0], di[1], di[2], di[3], di[4]))
			f.write('</dihedral>\n')
		f.write('</configuration>\n')
	def register(self, timestep):
		return
		