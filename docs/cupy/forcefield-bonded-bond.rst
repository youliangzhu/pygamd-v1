Bond interactions
-----------------

**Overview**

Bonds impose connected forces on specific pairs of particles to model chemical bonds.
The bonds are specified in :ref:`mst-format` configuration file with the format::

   bond
   bond_type(str)  particle_i(int)  particle_j(int)
   ...
   
By themselves, bonds do nothing. Only when a bond force object is instantiated in script(i.e. :py:class:`force.bond`), are bond forces actually calculated.

====================   ======================
:ref:`bond-function`   :py:class:`force.bond`
====================   ======================

.. image:: bond.png
    :width: 250 px
    :align: center
    :alt: Principle of bond stretching

.. _bond-function:
	
Bond functions
^^^^^^^^^^^^^^

Description:

   Function of bond interactions could be either the one called from bond interaction function libary, or the one defined by user himself.
   Bond interaction function libary contains harmonic function named as 'harmonic'.

   Harmonic function (harmonic)
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{bond}}(r) = \frac{1}{2}k\left( r - r_{0} \right)^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - spring constant ``k`` (in units of energy/distance^2)
    - :math:`r_0` - equilibrium length ``r0`` (in distance units)

.. py:class:: force.bond(info, func)

   Constructor of a bond interaction object.
 
   :param info: system information.
   :param func: function that is either a string or a device function. 
   
   .. py:function:: setParams(bond_type, param)
   
      specifies the bond interaction parameters with bond type and a list of parameters.

   Example::
   
      fb = pygamd.force.bond(info=mst, func='harmonic')
      fb.setParams(bond_type = 'A-A', param=[4.0, 0.0])#(param=[k, r0])
      fb.setParams(bond_type = 'A-B', param=[4.0, 0.0])#(param=[k, r0])
      fb.setParams(bond_type = 'B-B', param=[4.0, 0.0])#(param=[k, r0])
      app.add(fb)

Self-defined bond functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

   The device function for bond interactions could be written in script and conveyed 
   to kernel funciton for calculation.
   
   With the potential form of bond interactions :math:`p(r)`, the expression of parameters in script are: 

   * p = :math:`p(r)`
   * f = :math:`-(\partial p(r)/\partial r)(1/r)`  
   

   Function code template::

		@cuda.jit(device=True)
		def func(rsq, param, fp):
			p0 = param[0]
			p1 = param[1]
			...
			calculation codes
			...
			fp[0]=f
			fp[1]=p

		fb = pygamd.force.bond(info, func)
		fb.setParams(bond_type, param=[p0, p1, ...])
		app.add(fb)			
   
   Example::
   
		from numba import cuda
		import numba as nb

		@cuda.jit(device=True)
		def harmonic(rsq, param, fp):
			k = param[0]
			r0 = param[1]
			r = math.sqrt(rsq)
			f = k * (r0/r - nb.float32(1.0))
			p = nb.float32(0.5) * k * (r0 - r) * (r0 - r)
			fp[0]=f
			fp[1]=p

		fb = pygamd.force.bond(info=mst, func=harmonic)
		fb.setParams(bond_type='a-a', param=[100.0, 1.0])
		app.add(fb)	
 


	  
