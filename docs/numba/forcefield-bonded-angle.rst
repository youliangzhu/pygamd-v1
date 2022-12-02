Angle bending
-------------
**Overview**

Angles impose forces on specific triplets of particles to model chemical angles between two bonds.
The angles are specified in :ref:`mst-format` configuration file with the format::

   angle
   angle_type(str)  particle_i(int)  particle_j(int)  particle_k(int)
   ...
   
By themselves, angles do nothing. Only when an angle force object is instantiated(i.e. :py:class:`force.angle`), are angle forces actually calculated.

=====================   =======================
:ref:`angle-function`   :py:class:`force.angle`
=====================   =======================

.. image:: angle.png
    :width: 250 px
    :align: center
    :alt: Principle of angle bending

.. _angle-function:	
	
Angle functions
^^^^^^^^^^^^^^^

Description:

   Function of angle interactions could be either the one called from angle interaction function libary, or the one defined by user himself.
   Angle interaction function libary contains harmonic function named as 'harmonic' and harmonic cosine function named as 'harmonic_cos'.
   
   Harmonic function (harmonic)
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta) = \frac{1}{2}k\left( \theta -\theta_{0} \right)^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)
    - :math:`\theta_{0}` - equilibrium angle ``theta0`` (in radians)

    .. note::
	    Angles for the functions in library should be given in script in the unit of degree, and the program will convert them into radian automatically.

   Harmonic cosine function (harmonic_cos)		
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta)=k\left[ 1-\cos \left( \theta - {\theta}_{0} \right) \right]		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy)
    - :math:`\theta_{0}` - equilibrium angle ``theta0`` (in radians)
	
    .. note::
	    Angles for the functions in library should be given in script in the unit of degree, and the program will convert them into radian automatically.		

.. py:class:: force.angle(info, func)

   Constructor of an angle interaction object.
 
   :param info: system information.
   :param func: function that is either a string or a device function. 

   .. py:function:: setParams(angle_type, param)
   
      specifies the angle interaction parameters with angle type and a list of parameters.
	  
   Example::
   
      fa = pygamd.force.angle(info=mst, func='harmonic')
      fa.setParams(angle_type='a-a-a', param=[100.0, 90.0])
      app.add(fa)

Self-defined bond functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

   The device function for angle interactions could be written in script and conveyed 
   to kernel funciton for calculation.
   
   With the potential form of angle interactions :math:`p(\theta)`, the expression of parameters in script are: 

   * p = :math:`p(\theta)`
   * f = :math:`\partial p(\theta)/\partial \theta`  

   Function code template::

		@cuda.jit(device=True)
		def func(cos_abc, sin_abc, param, fp):
			p0 = param[0]
			p1 = param[1]
			...
			calculation codes
			...
			fp[0]=f
			fp[1]=p

		fa = pygamd.force.angle(info, func)
		fa.setParams(bond_type, param=[p0, p1, ...])
		app.add(fa)			
   
   Example::
   
		from numba import cuda
		import numba as nb

		@cuda.jit(device=True)
		def harmonic(cos_abc, sin_abc, param, fp):
			k = param[0]
			t0 = param[1]
			dth = math.acos(cos_abc) - math.pi*t0/180.0
			f = k * dth
			p = nb.float32(0.5) * f * dth
			fp[0]=f
			fp[1]=p
		
		fa = pygamd.force.angle(info=mst, func=harmonic)
		fa.setParams(angle_type='a-a-a', param=[400.0, 90.0])#param=[k, t0]
		app.add(fa)
  
	  