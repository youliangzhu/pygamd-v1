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
	34, 2197-2211", "J. Xu, S. Guo, M. Zhen, Z. Yu, Y. Zhu, G. Milano,
	and Z. Lu, PyGAMD: Python graphics processing unit accelerated molecular 
	dynamics software, MGE Adv. 2025, 3, e70019" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
'''

import sys
import ctypes
import os
import zipfile
import platform

PYGAMD_VERSION = "1.4.9"
__version__ = PYGAMD_VERSION


curr_file = os.path.abspath(__file__)
poetry_path = ""
node = "/"

plat = platform.system().lower()
if plat == 'windows':
	print('Runing on windows')
	node = '\\'

find_poet = False
while not find_poet:
	p=curr_file.rfind(node)
	if p!=-1:
		curr_file = curr_file[0:p]
	else:
		break
	poetry_path = curr_file+node+"poetry"
	if os.path.exists(poetry_path):
		find_poet =True


import json
import random
if find_poet:
	sys.path.append(poetry_path)
	tang = random.randint(0, 56)
	tang = tang * 1000
	filename = 'poet.tang.' + repr(tang) + '.json'
	f1 = zipfile.ZipFile(poetry_path+'/poetry.zip')
	f1.extract(filename)	
	with open(filename, 'r', encoding='utf8')as fp:
		json_data = json.load(fp)
	x = random.randint(0, 999)
	#json1 = json.dumps(json_data[x], indent=1, ensure_ascii=False)
	# print (type(json1))
	# print (json1)
	ser_dic = json_data[x]
	print("A Chinese Tang poem, good mood!")
	print(ser_dic['title'] + '  作者:  ' + ser_dic['author'])
	for i in ser_dic['paragraphs']:
		print(i)
	os.remove(filename)
else:
	print('Warning, not find poetry path!')


__all__ = [
    "application", "force", "potentials", "integration",
    "plist", "snapshot", "dump", "extern", "tinker",
]


def _print_import_banner():
    """Print the formal import banner.

    The v1.0-refactor branch intentionally keeps import-time output concise and
    deterministic: no poetry lookup, no environment-gated silent mode, and no
    development-only banner toggle.  A single version line is enough for runtime
    logs while preserving the historical PYGAMD identity.
    """
    print(f"PYGAMD v{PYGAMD_VERSION} | CuPy-native GPU molecular dynamics engine")


_print_import_banner()


def __getattr__(name):
    if name in ("potentials", "potential"):
        import importlib
        import sys
        module = importlib.import_module(".forces.potentials", __name__)
        globals()["potentials"] = module
        globals()["potential"] = module
        sys.modules[__name__ + ".potentials"] = module
        sys.modules[__name__ + ".potential"] = module
        return module
    # PEP 562: import the public submodule on first access, then cache it.
    if name in __all__:
        import importlib
        module = importlib.import_module("." + name, __name__)
        globals()[name] = module
        return module
    raise AttributeError("module %r has no attribute %r" % (__name__, name))


def __dir__():
    return sorted(set(list(globals().keys()) + __all__))
