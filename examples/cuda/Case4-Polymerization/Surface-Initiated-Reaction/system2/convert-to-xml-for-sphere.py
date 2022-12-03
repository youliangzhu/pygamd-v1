import os
import string
import random
#from __future__ import print_function
lx=30.0	     #box size
ly=30.0	
lz=30.0	
sigma=0.2    # surfacial density of initiators
radius=3     #radius of sphere
rmin=1.0     # the minimum distance between any two initiators

dat=file("icover.3.312.5.1.txt")   # downloaded from http://neilsloane.com/icosahedral.codes/index.html. 
data=[]
for line in dat:
		lin=line.strip('\n')
		data.append(float(lin))
nline=len(data)/3
pos=[]
type=[]
for j in range(0, nline):
	pos.append([data[j*3]*radius, data[j*3+1]*radius, data[j*3+2]*radius])
	type.append("A")

# distance check between initiators
npp=len(pos)
id_list=[]
def check(id):
	for i in range(0, len(id_list)):
		k=id_list[i]
		dx = pos[k][0]-pos[id][0]
		dy = pos[k][1]-pos[id][1]	
		dz = pos[k][2]-pos[id][2]
		rsq=dx*dx + dy*dy + dz*dz
		if(rsq<rmin*rmin):
			return False
	id_list.append(id)
	return True

for i in range(0, int(sigma*float(npp))):
	id=int(random.random()*float(npp))
	while (not check(id)):	
		id = int(random.random()*float(npp))
	type[id]="B"

# out put xml file
opd=open("sphere.xml","w")		
opd.write('<?xml version ="1.0" encoding ="UTF-8" ?>\n')
opd.write('<galamost_xml version="1.3">\n')
opd.write('<configuration time_step="0" dimensions="3" natoms="'+str(len(pos))+'" >\n')
opd.write('<box lx="'+str(lx)+'" ly="'+str(ly)+'" lz="'+str(lz)+'"/>\n')
opd.write('<position num="'+str(len(pos))+'">\n')
for i in range(0, len(pos)):
	opd.write(str(pos[i][0])+"  "+str(pos[i][1])+"  "+str(pos[i][2])+"\n")
opd.write('</position>\n')
opd.write('<type num="'+str(len(type))+'">\n')
for i in range(0, len(type)):
	opd.write(type[i]+"\n")
opd.write('</type>\n')
opd.write('</configuration>\n')
opd.write('</galamost_xml>\n')	
