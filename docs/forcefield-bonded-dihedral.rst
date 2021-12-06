Dihedral torsion
----------------
**Overview**

Dihedrals impose forces on specific quadruplets of particles to model the rotation about chemical bonds.
The dihedrals are specified in :ref:`mst-format` configuration file with the format::

   dihedral
   dihedral_type(str)  particle_i(int)  particle_j(int)  particle_k(int)  particle_l(int)
   ...
   
By themselves, dihedrals do nothing. Only when a dihedral force object is instantiated(i.e. :py:class:`force.dihedral`), are dihedral forces actually calculated.

========================   ==========================
:ref:`dihedral-function`   :py:class:`force.dihedral`
========================   ==========================

.. image:: dihedral.png
    :width: 250 px
    :align: center
    :alt: Principle of dihedral torsion

.. _dihedral-function:	

Dihedral functions
^^^^^^^^^^^^^^^^^^

Description:

   Function of angle interactions could be either the one called from angle interaction function libary, or the one defined by user himself.
   Angle interaction function libary contains harmonic function named as 'harmonic' and harmonic cosine function named as 'harmonic_cos'.

   Harmonic function for proper dihedrals(harmonic)
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)=k\left[ 1+f\cos \left( \varphi-\delta \right) \right]		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - multiplicative constant ``k`` (in units of energy)
    - :math:`\delta` - phase shift angle ``delta`` (in radians)
    - :math:`f` - factor ``f`` (unitless)	
      - *optional*: defaults to `-1.0`		
	
    .. note::
	    Dihedral angles for the functions in library should be given in script in the unit of degree, and the program will convert them into radian automatically.
		
   Harmonic function for improper dihedrals (harmonic)		
    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)=k\left( \varphi-\delta \right)^2		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)
    - :math:`\delta` - phase shift angle ``delta`` (in radians)

    .. note::
	    Dihedral angles for the functions in library should be given in script in the unit of degree, and the program will convert them into radian automatically.			

.. py:class:: force.dihedral(info, func)

   Constructor of a dihedral interaction object.
 
   :param info: system information.
   :param func: function that is either a string or a device function. 

   .. py:function:: setParams(dihedral_type, param, term='proper')
   
      specifies the dihedral interaction parameters with dihedral type, a list of parameters and the term of dihedral.
      The term candidates of dihedral are 'proper' and 'improper' with the default 'proper'. 	

   .. py:function:: setCosFactor(factor)
   
      specifies the factor of harmonic function for proper dihedral.		  
	  
   Example::
   
      fd = pygamd.force.dihedral(info=mst, func='harmonic')
      fd.setParams(dihedral_type='a-a-a-a', param=[100.0, 90.0])
      app.add(fd)	

Self-defined bond functions
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

   The device function for dihedral interactions could be written in script and conveyed 
   to kernel funciton for calculation.
   
   With the potential form of dihedral interactions :math:`p(\varphi)`, the expression of parameters in script are: 

   * p = :math:`p(\varphi)`
   * f = :math:`\triangle p(\varphi)/\triangle \varphi`     

   Function code template::

		@cuda.jit(device=True)
		def func(cos_abcd, sin_abcd, param, fp):
			p0 = param[0]
			p1 = param[1]
			...
			calculation codes
			...
			fp[0]=f
			fp[1]=p

		fd = pygamd.force.dihedral(info, func)
		fd.setParams(dihedral_type, param=[p0, p1, ...])
		app.add(fd)			
   
   Example::
   
		from numba import cuda
		import numba as nb
		
		@cuda.jit(device=True)
		def harmonic(cos_abcd, sin_abcd, param, fp):
			k = param[0]
			cos_phi0 = param[1]
			sin_phi0 = param[2]
			cos_factor = param[3]
			f = cos_factor * (-sin_abcd*cos_phi0 + cos_abcd*sin_phi0)
			p = nb.float32(1.0) + cos_factor * (cos_abcd*cos_phi0 + sin_abcd*sin_phi0)
			fp[0]=-k*f
			fp[1]=k*p

		fd = pygamd.force.dihedral(info=mst, func=harmonic)
		fd.setParams(dihedral_type='a-a-a-a', param=[100.0, math.cos(math.pi), math.sin(math.pi), -1.0])
		app.add(fd)	


