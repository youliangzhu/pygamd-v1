import jugamd
from numba import cuda
import numba as nb
import math


mst = jugamd.snapshot.read("example.mst")
app = jugamd.application.dynamics(info=mst, dt=0.001)

@cuda.jit(device=True)
def lj(rsq, param, fp):
	epsilon = param[0]
	sigma = param[1]
	alpha = param[2]
	rcut = param[3]
	if rsq<rcut*rcut:
		sigma2 = sigma*sigma
		r2inv = sigma2/rsq;
		r6inv = r2inv * r2inv * r2inv;
		f = nb.float32(4.0) * epsilon * r2inv * r6inv * (nb.float32(12.0) * r6inv - nb.float32(6.0) * alpha)/sigma2	
		p = nb.float32(4.0) * epsilon * r6inv * ( r6inv - nb.float32(1.0))
		fp[0]=f
		fp[1]=p	

fn = jugamd.force.nonbonded(info=mst, rcut=3.0, func=lj)
fn.setParams(type_i="A", type_j="A", param=[1.0, 1.0, 1.0, 3.0])
app.add(fn)

@cuda.jit(device=True)
def bond_harmonic(rsq, param, fp):
	k = param[0]
	r0 = param[1]
	r = math.sqrt(rsq)
	f = k * (r0/r - nb.float32(1.0))
	p = nb.float32(0.5) * k * (r0 - r) * (r0 - r)
	fp[0]=f
	fp[1]=p

fb = jugamd.force.bond(info=mst, func=bond_harmonic)
fb.setParams(bond_type = 'A-A', param=[1000.0, 0.5])#param=[k, r0]
app.add(fb)


@cuda.jit(device=True)
def angle_harmonic(cos_abc, sin_abc, param, fp):
	k = param[0]
	t0 = param[1]
	dth = math.acos(cos_abc) - math.pi*t0/nb.float32(180.0)
	f = k * dth
	p = nb.float32(0.5) * f * dth
	fp[0]=f
	fp[1]=p

fa = jugamd.force.angle(info=mst, func=angle_harmonic)
fa.setParams(angle_type='A-A-A', param=[100.0, 120.0])#param=[k, t0]
app.add(fa)

@cuda.jit(device=True)
def dihedral_harmonic(cos_abcd, sin_abcd, param, fp):
	k = param[0]
	cos_phi0 = param[1]
	sin_phi0 = param[2]
	cos_factor = param[3]
	f = cos_factor * (-sin_abcd*cos_phi0 + cos_abcd*sin_phi0)
	p = nb.float32(1.0) + cos_factor * (cos_abcd*cos_phi0 + sin_abcd*sin_phi0)
	fp[0]=-k*f
	fp[1]=k*p
	
fd = jugamd.force.dihedral(info=mst, func=dihedral_harmonic)
fd.setParams(dihedral_type = 'A-A-A-A', param=[25.0, math.cos(math.pi), math.sin(math.pi), -1.0])#param=[k, cos_phi0, sin_phi0, cos_factor]
app.add(fd)

nvt = jugamd.integration.nvt(info=mst, group='all', method="nh", tau=0.5, temperature=1.0)
app.add(nvt)

dd = jugamd.dump.data(info=mst, group='all', file='data.log', period=100)
app.add(dd)

dm = jugamd.dump.mst(info=mst, group='all', file='p.mst', period=10000)
app.add(dm)

#run 
app.run(20000)

