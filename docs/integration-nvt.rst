NVT ensemble
============

**Overview**

====================   ===========================
:ref:`nvt`             :py:class:`integration.nvt`
====================   ===========================


.. _nvt:

Thermostat
----------

.. py:class:: integration.nvt(info, group, method, tau, temperature)

   Constructor of a NVT NoseHoover thermostat object for a group of particles.
	  
   :param info: system information.
   :param group: a group of particles.	
   :param method: thermostat method, the candidates are "nh" for nose hoover.	   
   :param tau: thermostat coupling parameter.
   :param temperature: temperature.		  

   .. py:function:: setT(float T)
   
      specifies the temperature as a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature as a function of time steps.
	  
   Example::
   
      inn = pygamd.integration.nvt(info=mst, group=['a'], method="nh", tau=1.0, temperature=1.0)
      app.add(inn)
