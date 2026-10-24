NVT ensemble
============

**Overview**

====================   ============================
:ref:`nvt`             :py:class:`integration.nvt`
:ref:`gwvv`            :py:class:`integration.gwvv`
:ref:`bd`              :py:class:`integration.bd`
====================   ============================


.. _nvt:

NVT
---

.. py:class:: integration.nvt(info, group, method, tau, temperature)

   Constructor of a NVT thermostat object for a group of particles.
	  
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


.. _gwvv:

GWVV
----

.. py:class:: integration.gwvv(info, group)

   Constructor of a GWVV thermostat object for a group of particles.
	  
   :param info: system information.
   :param group: a group of particles.

   Example::
   
      inn = pygamd.integration.gwvv(info=mst, group='all')
      app.add(inn)

.. _bd:

BD
--

.. py:class:: integration.bd(info, group, temperature)

   Constructor of a Brownian Dynamics thermostat object for a group of particles.
	  
   :param info: system information.
   :param group: a group of particles.
   :param temperature: temperature.

   .. py:function:: setParams(string typ, float gamma)
   
      specifies the gamma parameter for particle type.
	  
   Example::
   
      inn = pygamd.integration.bd(info=mst, group="all", temperature=1.0)
      app.add(inn)

