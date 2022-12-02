NVT ensemble
============

**Overview**

====================   =========================
:ref:`nh-nvt`          :py:class:`NoseHooverNvt`
:ref:`berendsen-nvt`   :py:class:`BerendsenNvt`
:ref:`andersen-nvt`    :py:class:`AndersenNvt`
:ref:`bd-nvt`          :py:class:`BdNvt`
:ref:`nvt-rigid`       :py:class:`NvtRigid`
:ref:`bd-nvt-rigid`    :py:class:`BdNvtRigid`
====================   =========================


.. _nh-nvt:

Nose Hoover thermostat
----------------------

.. py:class:: NoseHooverNVT(all_info, group, comp_info, T, tauT)

   The constructor of a NVT NoseHoover thermostat object for a group of particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param ComputeInfo comp_info: The object of calculation of collective information.	   
   :param float T: The temperature.  
   :param float tauT: The thermostat coupling.		  

   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a defined varying value by time step.
	  
   Example::
   
      group = gala.ParticleSet(all_info, 'all')
      comp_info = gala.ComputeInfo(all_info, group)
      nh = gala.NoseHooverNVT(all_info, group, comp_info, 1.0, 0.5)
      app.add(nh)

.. _berendsen-nvt:

Berendsen thermostat
--------------------

.. py:class:: BerendsenNVT(all_info, group, comp_info, T, tauT)

   The constructor of a NVT Berendsen thermostat object for a group of particles.
	 
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param ComputeInfo comp_info: The object of calculation of collective information.	   
   :param float T: The temperature.  
   :param float tauT: The thermostat coupling parameter.	

   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
      
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time steps.
   
.. _andersen-nvt:
   
Andersen thermostat
-------------------

.. py:class:: AndersenNVT(all_info, group, T, gamma, seed)

   The constructor of a NVT Andersen thermostat object for a group of particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param float T: The temperature.  
   :param float gamma: The collision frequency.		  
   :param int seed: The seed of random number generator.	

   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time steps.
	  
   Example::
   
      an = gala.AndersenNVT(all_info,group,1.0,10.0, 12345)
      app.add(an)

.. _bd-nvt:	  
	  
Brownian dynamic thermostat
---------------------------

Description:

    The particles are integrated forward in time according to the Langevin equations of motion:

    .. math::

        m \frac{d\vec{v}}{dt} = \vec{F}_\mathrm{C} - \gamma \cdot \vec{v} + \vec{F}_\mathrm{R}

        \langle \vec{F}_\mathrm{R} \rangle = 0

        \langle |\vec{F}_\mathrm{R}|^2 \rangle = 2 d kT \gamma / \delta t
		
    - :math:`\gamma` - *gamma* (unitless) - *optional*: defaults to 1.0
	
    where :math:`\vec{F}_\mathrm{C}` is the force on the particle from all potentials and constraint forces,
    :math:`\gamma` is the drag coefficient, :math:`\vec{v}` is the particle's velocity, :math:`\vec{F}_\mathrm{R}`
    is a uniform random force, and :math:`d` is the dimensionality of the system (2 or 3).  The magnitude of
    the random force is chosen via the fluctuation-dissipation theorem to be consistent with the specified drag and temperature, :math:`T`.
    When :math:`kT=0`, the random force :math:`\vec{F}_\mathrm{R}=0`.

.. py:class:: BDNVT(all_info, group, T, seed)

   The constructor of a Brownian NVT thermostat object for a group of particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param float T: The temperature.    
   :param int seed: The seed of random number generator.		  

   .. py:function:: setGamma(float gamma)
   
      specifies the gamma with a constant value.
	  
   .. py:function:: setGamma(string type, float gamma)
   
      specifies the gamma of a particle type.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
 
      specifies the temperature with a varying value by time step.
	  
   Example::
   
      group = gala.ParticleSet(all_info, 'all')
      bdnvt = gala.BDNVT(all_info, group, 1.0, 123)
      app.add(bdnvt)

.. _nvt-rigid:

NVT for rigid body
------------------

.. py:class:: NVTRigid(AllInfo all_info, ParticleSet group, float T, float tauT)

   The constructor of a NVT thermostat object for rigid bodies.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param float T: The temperature.    
   :param float tauT: The thermostat coupling parameter.  

   .. py:function:: setT(float T)
   
      specifies the temperature with a fixed value.
	  
   .. py:function:: setT(Variant vT)
   
      pecifies the temperature with a varying value by time step.
	  
   Example::
   
      bgroup = gala.ParticleSet(all_info, 'body')
      rigidnvt = gala.NVTRigid(all_info, bgroup, 1.0, 10.0)
      app.add(rigidnvt)

.. _bd-nvt-rigid:	  
	  
Brownian dynamic for rigid body
-------------------------------

Please see :ref:`bd-nvt` for the theory.

.. py:class:: BDNVTRigid(all_info, group, T, seed)

   The constructor of a Brownian NVT thermostat object for rigid bodies.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param float T: The temperature.    
   :param int seed: The seed of random number generator.		  

   .. py:function:: setGamma(float gamma)
   
      specifies the gamma of Brownian method with a constant value.
	  
   .. py:function:: setGamma(const std::string & type, float gamma)
   
      specifies the gamma of Brownian method of a particle type.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time step.
	  
   Example::
   
      bgroup = gala.ParticleSet(all_info, 'body')
      bdrigidnvt = gala.BDNVTRigid(all_info, bgroup, 1.0, 123)
      app.add(bdrigidnvt)
