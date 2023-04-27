PYGAMD: Python GPU-Accelerated Molecular Dynamics Software
==========================================================

Version 1
---------

Copyright: Lu group
-------------------

Molecular dynamics is an essential tool in the scientific research on
colloids, polymers, biomolecules, and etc. ``PYGAMD`` is a platform
based on Python3 Numba where users can define their models and methods
by themselves. The simulation models could be ranging from all-atom to
coarse-grained scales.

Installation
------------

pygamd:

::

   python3 setup.py install

::

   pip install pygamd

Requirements:
~~~~~~~~~~~~~

1. Python3 including numba, numpy, and pybind11 packages
2. NVIDIA CUDA Toolkit >= 7.0

Citation
--------

::

   To cite PYGAMD in publications use:
    
   Thanks to Zhongyuan Lu group for providing PYGAMD software in our research work.

Documentation
-------------

Official website: `https://pygamd.com <https://pygamd.com>`__ Online
manual could be read here:
`https://pygamd-v1.readthedocs.io/en/latest/ <https://pygamd-v1.readthedocs.io/en/latest/>`__.
Tutorials written by jupyter notebook are given here:
`https://nbviewer.jupyter.org/github/youliangzhu/pygamd-v1/tree/main/tutorials/index.ipynb <https://nbviewer.jupyter.org/github/youliangzhu/pygamd-v1/tree/main/tutorials/index.ipynb>`__.
More examples could be found here:
`https://github.com/youliangzhu/pygamd-v1/tree/main/examples <https://github.com/youliangzhu/pygamd-v1/tree/main/examples>`__.

Example: DPD simulation of diblock copolymer
--------------------------------------------

1 First step: generate configuration

::

   from poetry import molgen

   mol1=molgen.Molecule(10)#particle number
   mol1.setParticleTypes("A,A,A,A,A,B,B,B,B,B")#type
   mol1.setTopology("0-1,1-2,2-3,3-4,4-5,5-6,6-7,7-8,8-9")#topology
   mol1.setBondLength(0.75)#bond length
   mol1.setMass(1.0)#mass


   gen=molgen.Generators(20,20,20) # box size in x, y, and z direction
   gen.addMolecule(mol1,2400)#molecule, the number of molecules
   gen.outPutMST("A5B5") #file name

2 Second step: run simulation

::

   import pygamd
       
   mst = pygamd.snapshot.read("A5B5.mst")
   app = pygamd.application.dynamics(info=mst, dt=0.04)

   fn = pygamd.force.dpd(info=mst, rcut=1.0)
   fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
   fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
   fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
   app.add(fn)

   fb = pygamd.force.bond(info=mst, func='harmonic')
   fb.setParams(bond_type = 'A-A', param=[4.0, 0.0])# param=[k, r0]
   fb.setParams(bond_type = 'A-B', param=[4.0, 0.0])# param=[k, r0]
   fb.setParams(bond_type = 'B-B', param=[4.0, 0.0])# param=[k, r0]
   app.add(fb)

   inn = pygamd.integration.gwvv(info=mst, group='all')
   app.add(inn)

   dm = pygamd.dump.mst(info=mst, group=['A', 'B'], file='p.mst', period=10000)
   app.add(dm)

   di = pygamd.dump.data(info=mst, group='all', file='data.log', period=100)
   app.add(di)

   app.run(500000)

Contributing
------------

We welcome contributions to PYGAMD. Whether it is reporting a bug,
starting a discussion by asking a question, or proposing/requesting a
new feature, please go by creating a new issue here
(`https://github.com/youliangzhu/pygamd-v1/issues/ <https://github.com/youliangzhu/pygamd-v1/issues/>`__)
or writing an email to the author Dr. You-Liang Zhu (Email:
ylzhu@pygamd.com) so that we can talk about it. Please note that this
project is released with a Contributor Code of Conduct. By participating
in this project you agree to abide by its terms.
