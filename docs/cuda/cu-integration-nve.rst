NVE ensemble
============

**Overview**

========================   =====================
:ref:`nve`                 :py:class:`NVE`
:ref:`nve-rigid`           :py:class:`NVERigid`
:ref:`nve-rigid-tunable`   :py:class:`TranRigid`
========================   =====================

.. _nve:

NVE thermostat
--------------

.. py:class:: NVE(all_info, group)

   The constructor of a NVE thermostat object for a group of particles.

   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.

   .. py:function:: setZeroForce(bool switch) 
   
      switches the function of making all force to be zero (the default is False).
   
   Example::
   
      group = gala.ParticleSet(all_info,'all')
      thermo = gala.NVE(all_info,group)
      app.add(thermo)

.. _nve-rigid:	  
	  
NVE for rigid body
------------------

.. py:class:: NVERigid(all_info, group)

   The constructor of a NVE thermostat object for rigid bodies.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.
   
   Example::
   
      bgroup = gala.ParticleSet(all_info, 'body')
      rigidnve = gala.NVERigid(all_info, bgroup)
      app.add(rigidnve)
	
.. _nve-rigid-tunable:
	
NVE for rigid body with tunable freedoms
----------------------------------------

.. py:class:: TranRigid(all_info, group)

   The constructor of a NVE thermostat object for rigid bodies for defined freedoms.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	  

   .. py:function:: setTraDimension(bool x, bool y, bool z)
   
      switches the freedoms of translocation in x y z directions.
	  
   .. py:function:: setRotDimension(bool x, bool y, bool z)
   
      switches the freedoms of rotation in x y z directions.
	  
   Example::
   
      rigidnve = gala.TranRigid (all_info, bgroup)
      rigidnve.setTraDimension(True, True, True)
      rigidnve.setRotDimension(True, True, True)
      app.add(rigidnve)
	  
