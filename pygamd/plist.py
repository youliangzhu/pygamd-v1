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

# CuPy-native main line: pygamd.plist.neighbor -> CuPy nlist/cplist -> CuPy clist
# (pygamd/plists/{nlist,cplist,clist}.py, pure cupy.RawKernel).  No Numba
# clist/nlist here; the legacy Numba versions live in archive_legacy/.
#
# Status: the CuPy nlist supports tag-space exclusions for bond / angle /
# dihedral connectivity.  Exclusions are applied while building the neighbor
# list, so RawPair and RawDPD consume the same filtered list without adding
# branches to force kernels.
from pygamd import plists
import os


_DEFAULT_CLUSTER_PAIR_THRESHOLD = 500000


def _read_cluster_threshold():
	raw = os.environ.get("PYGAMD_CLUSTER_PAIR_THRESHOLD", "")
	if raw == "":
		return _DEFAULT_CLUSTER_PAIR_THRESHOLD
	try:
		value = int(str(raw).strip())
	except (TypeError, ValueError):
		return _DEFAULT_CLUSTER_PAIR_THRESHOLD
	return max(0, value)


def _select_pairlist_algorithm(info):
	'''Return "nlist" or "cluster_pair".

	Default policy: use cluster pair list only when npa > 50000.  Override with
	PYGAMD_PAIRLIST_MODE=nlist|cluster|cluster_pair|auto for diagnostics.
	'''
	mode = os.environ.get("PYGAMD_PAIRLIST_MODE", "auto").strip().lower()
	if mode in ("nlist", "neighbor", "neighbor_list"):
		return "nlist"
	if mode in ("cluster", "cluster_pair", "cplist", "clusterpair"):
		return "cluster_pair"
	threshold = _read_cluster_threshold()
	return "cluster_pair" if int(info.npa) > threshold else "nlist"


def normalize_exclusion(exclusion):
	'''Normalize user exclusion spec to a deterministic tuple.

	Accepted forms:
	  None / [] / False       -> None
	  "bond" / "angle"     -> ("bond",) / ("angle",)
	  "bond,angle"          -> ("bond", "angle")
	  ["bond", "angle"]    -> ("bond", "angle")
	  "1-2" / "1-3" / "1-4" -> bond / angle / dihedral
	  "all" / "bonded"     -> bond + angle + dihedral

	Angle exclusions mean the 1-3 endpoint pair of an angle.  Dihedral exclusions
	mean the 1-4 endpoint pair.  Bond pairs are not implied by angle/dihedral; pass
	["bond", "angle"] when both are desired.
	'''
	if exclusion is None or exclusion is False:
		return None
	if isinstance(exclusion, str):
		items = []
		for tok in exclusion.replace(';', ',').replace('+', ',').split(','):
			items.extend(tok.strip().split())
	else:
		try:
			items = list(exclusion)
		except TypeError:
			items = [exclusion]
	aliases = {
		"bond": "bond", "bonds": "bond", "1-2": "bond", "12": "bond",
		"angle": "angle", "angles": "angle", "1-3": "angle", "13": "angle",
		"dihedral": "dihedral", "dihedrals": "dihedral", "1-4": "dihedral", "14": "dihedral",
	}
	order = ["bond", "angle", "dihedral"]
	seen = set()
	for item in items:
		if item is None or item is False:
			continue
		name = str(item).strip().lower().replace('_', '-').replace(' ', '')
		if name in ("", "none", "false", "off", "0"):
			continue
		if name in ("all", "bonded", "full"):
			seen.update(order)
			continue
		if name not in aliases:
			raise RuntimeError("Unsupported exclusion kind %r; supported: bond, angle, dihedral, all" % (item,))
		seen.add(aliases[name])
	result = tuple(k for k in order if k in seen)
	return result if result else None


class neighbor:
	def __init__(self, info, rcut, exclusion):
		self.rcut=rcut
		self.exclusion = normalize_exclusion(exclusion)
		# rbuff = 0.1 * rcut (unchanged).
		self.algorithm = _select_pairlist_algorithm(info)
		if self.algorithm == "cluster_pair":
			self.data = plists.cplist.cplist(info, rcut, rcut*0.1, self.exclusion)
		else:
			self.data=plists.nlist.nlist(info, rcut, rcut*0.1, self.exclusion)
		self.last_ts = 0xffffffff

	def calculate(self, timestep):
		if timestep == self.last_ts:
			return
		self.last_ts = timestep
		self.data.calculate(False, timestep)

	def output(self):
		self.data.speak()
