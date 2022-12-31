import galamost
import math
	
	
def pair(width, func, rmin, rmax, coeff):
	ptable = galamost.vector_real2()
	dr = rmax/width
	for i in range(0, width):
		r = dr * i
		if r<rmin:
			potential = func(rmin, **coeff)
		else:
			potential = func(r, **coeff)
		ptable.append(galamost.ToReal2(r, potential))
	return ptable
	
def bond(width, func, rmin, rmax, coeff):
	ptable = galamost.vector_real2()
	dr= rmax/width
	for i in range(0, width):
		r = dr * i
		if r<rmin:
			potential = func(rmin, **coeff)
		else:
			potential = func(r, **coeff)
		ptable.append(galamost.ToReal2(r, potential))
	return ptable
	
def angle(width, func, coeff):
	ptable = galamost.vector_real2()
	dth = math.pi/width
	for i in range(0, width):
		th = dth * i
		potential = func(th, **coeff)
		ptable.append(galamost.ToReal2(th, potential))
	return ptable
	
def dihedral(width, func, coeff):
	ptable = galamost.vector_real2()
	dth = 2.0*math.pi/width
	for i in range(0, width):
		th = dth * i
		potential = func(th, **coeff)
		ptable.append(galamost.ToReal2(th, potential))
	return ptable	
	

