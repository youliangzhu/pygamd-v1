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

from pathlib import Path

from setuptools import setup, find_packages

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(name="pygamd",
    version="1.4.9",
    description="Python GPU-Accelerated Molecular Dynamics Software",
    long_description=long_description,
    long_description_content_type="text/x-rst",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3.12",
        "Intended Audience :: Science/Research",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    keywords="molecular dynamics",
    url="https://github.com/ylzhu/pygamd.py",
    author="You-Liang Zhu",
    author_email="ylzhu@pygamd.com",
    license="GPL-3.0-or-later",
    license_files=["LICENSE"],
    packages=find_packages(),
    python_requires=">=3.12,<3.13",
    install_requires=[
        "numpy>=1.26",
        "networkx>=3.0",
    ],
    extras_require={
        # CuPy wheels are CUDA-major-specific.  Keeping them as mutually
        # exclusive extras avoids a slow and usually unsuccessful source build
        # of the generic ``cupy`` package.
        "cuda12": ["cupy-cuda12x>=14.1,<15"],
        "cuda13": ["cupy-cuda13x>=14.1,<15"],
        "aa2cg": ["MDAnalysis>=2.9,<3"],
        "rdf": ["MDAnalysis>=2.9,<3"],
        "test": ["pytest>=8"],
    },
    entry_points={
        "console_scripts": [
            "multiscale=multiscale.cli:main",
        ],
    },
    package_data={
        "multiscale": ["py.typed"],
        "multiscale.forcefield.assets.opls": [
            "*.txt",
            "*.itp",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

