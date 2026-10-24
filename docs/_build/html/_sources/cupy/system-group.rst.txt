Group
=====
   
A group specifies the particles for certain functions, such as integration and data output. Group objects
will not be defined in script. Instead, they will be defined in other objects. And then, ususally only keywords or a list of
particle types and particle indexes are needed to indicate the particles.  
      
.. py:class:: chare.particle_set(info, group)

   :param info: system information.
   :param group: either a string or a python list. The string is a keyword with candidates "all", "body", "charge", and "nonbody". 
				The list could contain particle types and particle indexes.   

   Example::
   
      dm = pygamd.dump.mst(info=mst, group='all', file='p.mst', period=100000)
      app.add(dm)
	  
      dm = pygamd.dump.mst(info=mst, group=['A', 'B'], file='p.mst', period=100000)
      app.add(dm)


      dm = pygamd.dump.mst(info=mst, group=['A', 'B', 0, 1, 2], file='p.mst', period=100000)
      app.add(dm)	  
	  

   
   
   
   