NVT ensemble
============

**Overview**

=========================      ============================
:ref:`nh-nvt`                  :py:class:`NoseHooverNVT`
:ref:`nhc-nvt`                 :py:class:`NoseHooverChainNVT`
:ref:`berendsen-nvt`           :py:class:`BerendsenNVT`
:ref:`andersen-nvt`            :py:class:`AndersenNVT`
:ref:`langevin-nvt`            :py:class:`LangevinNVT`
:ref:`nvt-rigid`               :py:class:`NVTRigid`
:ref:`langevin-nvt-rigid`      :py:class:`LangevinNVTRigid`
=========================      ============================


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

.. _nhc-nvt:

Nose Hoover chain thermostat
----------------------------

.. py:class:: NoseHooverChainNVT(all_info, group, comp_info, T, tauT)

   The constructor of a NVT NoseHoover chain thermostat object for a group of particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param ComputeInfo comp_info: The object of calculation of collective information.	   
   :param float T: The temperature.  
   :param float tauT: The thermostat coupling.		  

   .. py:function:: setTau(Real tauT)
   
      specifies the thermostat coupling.
	  
   Example::
   
      group = gala.ParticleSet(all_info, 'all')
      comp_info = gala.ComputeInfo(all_info, group)
      nhc = gala.NoseHooverChainNVT(all_info, group, comp_info, 1.0, 0.5)
      app.add(nhc)


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

.. _langevin-nvt:	  
	  
Langevin dynamic thermostat
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

.. py:class:: LangevinNVT(all_info, group, T, seed)

   The constructor of a Langevin NVT thermostat object for a group of particles.
	  
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
      lnvt = gala.LangevinNVT(all_info, group, 1.0, 123)
      app.add(lnvt)

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

.. _langevin-nvt-rigid:	  
	  
Langevin dynamic for rigid body
-------------------------------

Please see :ref:`langevin-nvt` for the theory.

.. py:class:: LangevinNVTRigid(all_info, group, T, seed)

   The constructor of a Langevin NVT thermostat object for rigid bodies.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param float T: The temperature.    
   :param int seed: The seed of random number generator.		  

   .. py:function:: setGamma(float gamma)
   
      specifies the gamma of Langevin method with a constant value.
	  
   .. py:function:: setGamma(const std::string & type, float gamma)
   
      specifies the gamma of Langevin method of a particle type.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time step.
	  
   Example::
   
      bgroup = gala.ParticleSet(all_info, 'body')
      lrigidnvt = gala.LangevinNVTRigid(all_info, bgroup, 1.0, 123)
      app.add(lrigidnvt)
