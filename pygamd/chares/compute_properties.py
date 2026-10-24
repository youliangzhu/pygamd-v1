'''
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or
	modify it under the terms of the GNU General Public License.
	See the General Public License v3 for more details.
CORRESPONDENCE
	Dr. You-Liang Zhu  Email: ylzhu@pygamd.com
'''

# -----------------------------------------------------------------------------
# CuPy-native group property reductions (MVP).  No numba.  Consumes the CuInfo
# CuPy device arrays directly and returns the host-side sums dump.data needs:
#   momentum, temperature, pressure, potential.
#
# group_sums(cu, d_member) returns, for the group's member indices:
#   ke2_sum  = sum_i m_i * (vx^2 + vy^2 + vz^2)   (== 2 * kinetic energy)
#   vir_sum  = sum_i virial_potential[i, 0]
#   pot_sum  = sum_i virial_potential[i, 1]
#   mx,my,mz = sum_i m_i * v_i                     (group momentum components)
# computed with CuPy reductions; one device->host transfer per call.
# cupy is imported lazily; importing this module never requires a CUDA driver.
# -----------------------------------------------------------------------------

from pygamd.backend.cupy import import_cupy


def group_sums(cu, d_member):
	'''CuPy reduction over the group members.  `cu` is a CuInfo; `d_member` is a
	CuPy int32 array of member particle indices.  Returns a 6-tuple of python
	floats: (ke2_sum, vir_sum, pot_sum, mx, my, mz).'''
	cp = import_cupy("CuPy compute_properties")
	vm = cu.d_vel[d_member]                 # (nme, 4): vx, vy, vz, mass
	vpm = cu.d_virial_potential[d_member]   # (nme, 2): virial, potential
	vx = vm[:, 0]; vy = vm[:, 1]; vz = vm[:, 2]; m = vm[:, 3]
	ke2 = cp.sum(m * (vx * vx + vy * vy + vz * vz))
	vir = cp.sum(vpm[:, 0])
	pot = cp.sum(vpm[:, 1])
	mx = cp.sum(m * vx)
	my = cp.sum(m * vy)
	mz = cp.sum(m * vz)
	# single device->host transfer
	out = cp.asnumpy(cp.stack([ke2, vir, pot, mx, my, mz]))
	return (float(out[0]), float(out[1]), float(out[2]),
	        float(out[3]), float(out[4]), float(out[5]))
