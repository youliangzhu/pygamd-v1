Information computation
=======================

.. py:class:: ComputeInfo(all_info, group)

   The constructor of an object of computing some important information, including temperature, pressure, momentum, and potential of a group of particles.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.
   
   .. py:function:: setNdof(unsigned int nfreedom)
   
      sets the degree of freedom.
	  
   Example::
   
      comp_info = gala.ComputeInfo(all_info, group)

