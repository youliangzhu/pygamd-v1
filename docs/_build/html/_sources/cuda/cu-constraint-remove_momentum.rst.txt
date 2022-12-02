Remove CM momentum 
==================

.. py:class:: ZeroMomentum(all_info)

   The constructor of an object of removing the momentum of center mass of all particles.
   
   :param AllInfo all_info: The system information.
	   
   
.. py:class:: ZeroMomentum(all_info, group)

   specifies the method of removing the momentum of center mass of a group of particles.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of charged particles.    
   
   Example::
   
      zm = gala.ZeroMomentum(all_info)
      zm.setPeriod(10)
      app.add(zm)


