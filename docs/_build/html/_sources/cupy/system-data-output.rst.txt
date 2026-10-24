Data output
===========

Collective information
----------------------

.. py:class:: dump.data(info, group, file, period)

   Constructor of an information dump object for a group of particles.
   
   :param info: system information.	
   :param group: a group of particles.	   
   :param file: the name of output file.  
   :param period: the period of data output.
  
	  
   Example::
   
		dd = pygamd.dump.data(info=mst, group=['a'], file='data.log', period=100)
		app.add(dd)

MST dump
--------

.. py:class:: dump.mst(info, group, file, period, properties=None, split=False)

   Constructor of an object to dump MST files.

   :param info: system information.	
   :param group: a group of particles.	   
   :param file: the name of output file.  
   :param period: the period of data output.
   :param properties: the properties for output, candidates are
					'position', 'type', 'velocity', 'mass', 'image', 'force',
					'potential', 'virial', 'bond', 'angle', 'dihedral'
   :param split: if the trajectory file is splited into separated snapshot files    

   Example::
   
      dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=100000)
      app.add(dm)

XML dump
--------

.. py:class:: dump.xml(info, group, file, period, properties=None, split=False)

   Constructor of an object to dump XML files.

   :param info: system information.	
   :param group: a group of particles.	   
   :param file: the name of output file.  
   :param period: the period of data output.
   :param properties: the properties for output, candidates are
					'position', 'type', 'velocity', 'mass', 'image', 'force',
					'potential', 'virial', 'bond', 'angle', 'dihedral'
   :param split: if the trajectory file is splited into separated snapshot files    

   Example::
   
      dx = pygamd.dump.xml(info=mst, group='all', file='p', period=100000)
      app.add(dx)