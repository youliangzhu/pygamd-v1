Angle bending
-------------
**Overview**

Angles impose forces on specific triplets of particles to model chemical angles between two bonds.
The angles are specified in :ref:`xml-format` configuration file with the format::

   <angle>
   angle_type(str)  particle_i(int)  particle_j(int)  particle_k(int)
   ...
   </angle>
   
By themselves, angles do nothing. Only when you specify an angle force in script(i.e. :py:class:`AngleForceHarmonic`), are forces actually calculated between the listed particles.

============================   =================================
:ref:`harmonic-angle`          :py:class:`AngleForceHarmonic`
:ref:`harmonic-cosine-angle`   :py:class:`AngleForceHarmonicCos`
:ref:`cosine-angle`            :py:class:`AngleForceCos`
:ref:`LnExp-angle`             :py:class:`AngleForceLnExp`
============================   =================================

.. image:: angle.png
    :width: 250 px
    :align: center
    :alt: Principle of angle bending

.. _harmonic-angle:	
	
Harmonic angle potential
^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta) = \frac{1}{2}k\left( \theta -\theta_{0} \right)^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)
    - :math:`\theta_{0}` - equilibrium angle ``theta0`` (in radians)

    .. note::
	    The angles set in script are in the unit of degree, and the program will convert them into radian automatically.

.. py:class:: AngleForceHarmonic(all_info)

   The constructor of angle harmonic interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string type, float k, float theta0)
   
      specifies the angle harmonic force parameters with angle type, potential constant, and equilibrium angle degree.
	  
   Example::
   
      angleforce = gala.AngleForceHarmonic(all_info)
      angleforce.setParams('P-G-G', 25.000, 120.000)
      app.add(angleforce)

.. _harmonic-cosine-angle:
	  
Harmonic cosine angle potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta)=\frac{1}{2}k\left[\cos \left( \theta \right)-\cos \left( \theta_{0} \right)\right]^{2}	
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy)
    - :math:`\theta_{0}` - equilibrium angle ``theta0`` (in radians)
	
    .. note::
	    The angles set in script are in the unit of degree, and the program will convert them into radian automatically.

.. py:class:: AngleForceHarmonicCos(all_info)

   The constructor of angle cosine harmonic interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string type, float k, float theta0)
   
      specifies the angle cosine harmonic force parameters with angle type, potential constant, and equilibrium angle degree.
	  
   Example::
   
      angleforce = gala.AngleForceHarmonicCos(all_info)
      angleforce.setParams('P-G-G',25.000, 120.000)
      app.add(angleforce)

.. _cosine-angle:
	  
Cosine angle potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta)=k\left[ 1-\cos \left( \theta - {\theta}_{0} \right) \right]		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy)
    - :math:`\theta_{0}` - equilibrium angle ``theta0`` (in radians)
	
    .. note::
	    The angles set in script are in the unit of degree, and the program will convert them into radian automatically.	

.. py:class:: AngleForceCos(all_info)

   The constructor of angle cosine interaction object.
 
   :param AllInfo all_info: The system information.


   .. py:function:: setParams(string type, float k, float theta0)
   
      specifies the angle cosine force parameters with angle type, spring constant, and equilibrium angle degree.
	  
   Example::
   
      angleforce = gala.AngleForceCos(all_info)
      angleforce.setParams('P-G-G', 25.000, 120.000)
      app.add(angleforce)
  
 .. _LnExp-angle: 
  
LnExp angle potential
^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{angle}}(\theta)= -\frac{1}{2} k \log \left[ A \exp \left( -k_1 \left( \theta - {\theta}_{1} \right)^2 \right) + B \exp \left( -k_2 \left( \theta - {\theta}_{2} \right)^2 \right) \right]		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy) 
    - :math:`k_1, k_2, A, B` - potential parameters ``k1, k2, A, B``
    - :math:`\theta_{1}` - equilibrium angle ``theta1`` (in radians)
    - :math:`\theta_{2}` - equilibrium angle ``theta2`` (in radians)
	
    .. note::
	    The angles set in script are in the unit of degree, and the program will convert them into radian automatically.	

.. py:class:: AngleForceLnExp(all_info)

   The constructor of angle cosine interaction object.
 
   :param AllInfo all_info: The system information.


   .. py:function:: setParams(string type, float k, float k1, float k2, float theta1, float theta2, float A, float B)
   
      specifies the angle cosine force parameters with: angle type, spring constant, exponential factor1, exponential factor2,
      equilibrium angle degree1, equilibrium angle degree2, parameter A, and parameter B.
	  
   Example::
   
      angleforce = gala.AngleForceLnExp(all_info)
      angleforce.setParams('P-G-G', 25.000, 1.0, 1.0, 90.0, 180.0, 3.0, 2.0)
      app.add(angleforce)
  	  