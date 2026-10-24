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
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''

# -----------------------------------------------------------------------------
# CuPy-native MVP public force entry.
#
# The only supported pair-force path is the CuPy RawKernel backend
# (forces.pair.pair -> forces.pair).  The legacy Numba device-function pair
# path and the legacy bond/angle/dihedral/threebody/DPD/SLJ/DPK force classes
# have been removed from the main line (see archive_legacy/ for reference).
#
# User entry for raw non-charged pair potentials:
#     gem = pygamd.potentials.gem4_raw(fast_math=True)  # or lj_raw / pair_raw
#     fn  = pygamd.force.nonbonded(info=mst, rcut=2.5, func=gem)
# Built-in charged and bonded potentials use pygamd-numba-style string selectors:
#     pygamd.force.nonbonded_c(info=mst, rcut=1.2, func="lj_coulomb")
#     pygamd.force.nonbonded_c(info=mst, rcut=1.2, func="lj_coulomb_erfc")
#     pygamd.force.pppm(info=mst, rcut=1.2, kappa=2.5, epsilonr=15.0)
#     pygamd.force.bond(info=mst, func="harmonic")
#     pygamd.force.angle(info=mst, func="harmonic_cos")
#     fn.setParams(type_i="A", type_j="A", param=[...])
#     app.add(fn)
# `pygamd.force.pair` is an alias of `pygamd.force.nonbonded`.
# -----------------------------------------------------------------------------

from pygamd import forces
from pygamd import plist


class nonbonded:
	'''Non-bonded pair force (CuPy RawKernel).  `func` must be a RawPairPotential
	built with pygamd.potentials.pair_raw / lj_raw / gem4_raw.'''

	def __init__(self, info, rcut, func, exclusion=None):
		exclusion = plist.normalize_exclusion(exclusion)
		nl = info.find_plist(rcut, exclusion)
		if nl is None:
			nl = plist.neighbor(info, rcut, exclusion)
			info.plist.append(nl)
		self.info = info
		self.data = forces.pair.pair(info, nl, func)
		self.name = "force"
		self.subname = "nonbonded"
		self.check = [False for i in range(self.info.ntypes * self.info.ntypes)]
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
					raise RuntimeError('Error! the parameters between type ', self.info.typemap[i],
					                   ' and type ', self.info.typemap[j], ' has not been set!')


class nonbonded_c:
	'''Charged non-bonded pair force (CuPy RawKernel).

	Use func='lj_coulomb' for the pygamd-numba-compatible cutoff Coulomb,
	func='lj_coulomb_erfc' / 'coulomb_erfc' as the real-space partner of
	pygamd.force.pppm(...), or pass a custom pygamd.potentials.pair_charge(...).
	'''

	def __init__(self, info, rcut, func, exclusion=None):
		exclusion = plist.normalize_exclusion(exclusion)
		nl = info.find_plist(rcut, exclusion)
		if nl is None:
			nl = plist.neighbor(info, rcut, exclusion)
			info.plist.append(nl)
		self.info = info
		self.data = forces.pair.pair_charge(info, nl, func)
		self.name = "force"
		self.subname = "nonbonded_c"
		self.check = [False for i in range(self.info.ntypes * self.info.ntypes)]
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
					raise RuntimeError('Error! the parameters between type ', self.info.typemap[i],
					                   ' and type ', self.info.typemap[j], ' has not been set!')


# `pygamd.force.pair(...)` is the same CuPy RawKernel pair force as nonbonded.
pair = nonbonded


class pppm:
	'''PPPM reciprocal-space electrostatics (CuPy/cuFFT).

	Pair it with nonbonded_c(func='lj_coulomb_erfc' or 'coulomb_erfc') using the
	same kappa/epsilonr.  The force is additive and accumulates the long-range
	Coulomb contribution into the global force buffer.
	'''

	def __init__(self, info, grid=None, mesh_spacing=1.0, kappa=None, rcut=None,
	             epsilonr=1.0, coulomb_const=138.935, exclusion=None,
	             block_size=128, deconvolve=True):
		exclusion = plist.normalize_exclusion(exclusion)
		self.info = info
		self.data = forces.pppm.pppm(
			info=info,
			grid=grid,
			mesh_spacing=mesh_spacing,
			kappa=kappa,
			rcut=rcut,
			epsilonr=epsilonr,
			coulomb_const=coulomb_const,
			exclusion=exclusion,
			block_size=block_size,
			deconvolve=deconvolve,
		)
		self.name = "force"
		self.subname = "pppm"

	def compute(self, timestep):
		self.data.calculate(timestep)


# -----------------------------------------------------------------------------
# CuPy-native bonded forces.  These consume the CuInfo bond/angle/dihedral
# topology and accumulate into d_force / d_virial_potential (additive, after the
# global CuPy zero-force).  Built-in potentials are selected with strings, and
# scripts can pass custom CUDA device potentials created with
# pygamd.potentials.bond_raw / angle_raw / dihedral_raw.
#   bond:      built-in 'harmonic' or custom func(rsq, param, fp)
#   angle:     built-in 'harmonic'/'harmonic_cos' or custom func(cos, sin, param, fp)
#   dihedral:  built-in 'periodic'/'harmonic' or custom func(cos, sin, param, fp)
# -----------------------------------------------------------------------------

class bond:
	'''Bond force (CuPy RawKernel): built-in harmonic or custom device potential.'''

	def __init__(self, info, func="harmonic"):
		self.info = info
		self.data = forces.bond.bond(info, func)
		self.name = "force"
		self.subname = "bond"
		self.check = [False for i in range(self.data.ntypes)]
		self.first_compute = True

	def setParams(self, bond_type, param):
		self.data.setParams(bond_type, param, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)

	def check_parameters(self):
		tm = self.data.topo.typemap
		for tid in range(len(tm)):
			if not self.check[tid]:
				raise RuntimeError("Error! parameters for bond type '%s' have not been set!" % tm[tid])


class angle:
	'''Angle force (CuPy RawKernel): built-in harmonic variants or custom device potential.'''

	def __init__(self, info, func="harmonic"):
		self.info = info
		self.data = forces.angle.angle(info, func)
		self.name = "force"
		self.subname = "angle"
		self.check = [False for i in range(self.data.ntypes)]
		self.first_compute = True

	def setParams(self, angle_type, param):
		self.data.setParams(angle_type, param, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)

	def check_parameters(self):
		tm = self.data.topo.typemap
		for tid in range(len(tm)):
			if not self.check[tid]:
				raise RuntimeError("Error! parameters for angle type '%s' have not been set!" % tm[tid])


class dihedral:
	'''Dihedral force (CuPy RawKernel): built-in periodic/harmonic or custom device potential.'''

	def __init__(self, info, func="periodic"):
		self.info = info
		self.data = forces.dihedral.dihedral(info, func)
		self.name = "force"
		self.subname = "dihedral"
		self.check = [False for i in range(self.data.ntypes)]
		self.first_compute = True

	def setParams(self, dihedral_type, param, term="proper"):
		self.data.setParams(dihedral_type, param, term, self.check)

	def compute(self, timestep):
		if self.first_compute:
			self.check_parameters()
			self.first_compute = False
		self.data.calculate(timestep)

	def check_parameters(self):
		tm = self.data.topo.typemap
		for tid in range(len(tm)):
			if not self.check[tid]:
				raise RuntimeError("Error! parameters for dihedral type '%s' have not been set!" % tm[tid])


class dpd:
	'''DPD pairwise force (CuPy RawKernel): conservative + dissipative + random
	(Saru RNG).  setParams(type_i, type_j, alpha, sigma[, r_cut]).  Shares the
	CuPy neighbor list with RawPair.  Ordinary DPD uses velocity_source='vel';
	integration.dpdgwvv switches it to the predictor velocity_source='velo'.'''

	def __init__(self, info, rcut=1.0, temperature=1.0, seed=12345, exclusion=None):
		exclusion = plist.normalize_exclusion(exclusion)
		nl = info.find_plist(rcut, exclusion)
		if nl is None:
			nl = plist.neighbor(info, rcut, exclusion)
			info.plist.append(nl)
		self.info = info
		self.data = forces.dpd.dpd(info, nl, rcut, temperature, seed)
		self.name = "force"
		self.subname = "dpd"
		self.check = [False for i in range(self.info.ntypes * self.info.ntypes)]
		self.first_compute = True

	def setParams(self, type_i, type_j, alpha, sigma, r_cut=None):
		self.data.setParams(type_i, type_j, alpha, sigma, r_cut, self.check)

	def setT(self, temperature):
		self.data.setT(temperature)

	def setConsiderDiameter(self, con_dia, core_repulsion=0.0):
		self.data.setConsiderDiameter(con_dia, core_repulsion)

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
					raise RuntimeError("Error! DPD parameters between type '%s' and '%s' have not been set!"
					                   % (self.info.typemap[i], self.info.typemap[j]))
