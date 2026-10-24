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

# -----------------------------------------------------------------------------
# CuPy-native MVP public integration entry.
#
# Supported CuPy-native integrators are NVT, NVE, and BD, implemented through
# integrations.nvt / nve / bd.  Legacy Numba integrators are not part
# of the eager-graph primary route.
#
# User entries:
#     inn = pygamd.integration.nvt(info=mst, group="all", method="nh",
#                                  tau=0.5, temperature=0.02)
#     inn = pygamd.integration.nve(info=mst, group="all")
#     inn = pygamd.integration.bd(info=mst, group="all", temperature=0.02)
#     app.add(inn)
# -----------------------------------------------------------------------------

from pygamd import integrations
from pygamd import chare
from pygamd.backend.cupy.runtime import DPDGWVV_DEFAULT_LAMBDA


def _comp_info(info, group):
	ci = info.find_comp_info(group)
	if ci is None:
		ci = chare.comp_info(info, group)
		info.comp_info.append(ci)
	return ci


class nvt:
	'''Nose-Hoover NVT (CuPy RawKernel).'''

	def __init__(self, info, group, method, tau, temperature):
		if method != "nh":
			raise RuntimeError("Error: pygamd.integration.nvt only supports method='nh' "
			                   "(CuPy RawKernel Nose-Hoover) in the CuPy-native MVP.")
		self.data = integrations.nvt.nose_hoover(info, _comp_info(info, group), tau, temperature)
		self.name = "integration"
		self.subname = "nvt"

	def setT(self, temperature):
		self.data.setT(temperature)


class nve:
	'''Velocity-Verlet NVE (CuPy RawKernel).'''

	def __init__(self, info, group):
		self.data = integrations.nve.nve(info, _comp_info(info, group))
		self.name = "integration"
		self.subname = "nve"


class gwvv:
	"""Minimal Groot-Warren DPD velocity-Verlet integrator (CuPy RawKernel).

	This is the recommended public DPDGWVV entry.  Pass the DPD force wrapper so
	its RawDPD backend reads the predictor velocity field d_velo.  No separate legacy DPD integrator alias is exposed.
	"""

	def __init__(self, info, group, dpd_force=None, lambda_gw=DPDGWVV_DEFAULT_LAMBDA):
		self.data = integrations.dpd.gwvv(info, _comp_info(info, group), dpd_force=dpd_force, lambda_gw=lambda_gw)
		self.name = "integration"
		self.subname = "gwvv"

	def attachDPDForce(self, dpd_force):
		self.data.attach_dpd_force(dpd_force)

	def setLambda(self, lambda_gw):
		self.data.setLambda(lambda_gw)


class bd:
	'''Brownian / Langevin dynamics (CuPy RawKernel, Saru RNG).
	setParams(type, gamma) sets per-type friction; setT(value) the temperature.'''

	def __init__(self, info, group, temp, seed=12345):
		self.data = integrations.bd.bd(info, _comp_info(info, group), temp, seed)
		self.name = "integration"
		self.subname = "bd"

	def setParams(self, typ, gamma):
		self.data.setParams(typ, gamma)

	def setT(self, temperature):
		self.data.setT(temperature)
