Dihedral torsion
----------------
**Overview**

Dihedrals impose forces on specific quadruplets of particles to model the rotation about chemical bonds.
The dihedrals are specified in :ref:`xml-format` configuration file with the format::

   <dihedral>
   dihedral_type(str)  particle_i(int)  particle_j(int)  particle_k(int)  particle_l(int)
   ...
   </dihedral>
   
By themselves, dihedrals do nothing. Only when you specify a dihedral force in script(i.e. :py:class:`DihedralForceHarmonic`), are forces actually calculated between the listed particles.

=================================   ==========================================
:ref:`harmonic-dihedral`            :py:class:`DihedralForceHarmonic`
:ref:`harmonic-dihedral-improper`   :py:class:`DihedralForceHarmonic`
:ref:`opls-dihedral`                :py:class:`DihedralForceOplsCosine`
:ref:`rb-dihedral`                  :py:class:`DihedralForceRyckaertBellemans`
:ref:`amber-dihedral`               :py:class:`DihedralForceAmberCosine`
=================================   ==========================================

.. image:: dihedral.png
    :width: 250 px
    :align: center
    :alt: Principle of dihedral torsion

.. _harmonic-dihedral:	

Harmonic dihedral potential (proper)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

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
	    The dihedral angle delta in script are in the unit of degree, and the program will convert them into radian automatically.

.. py:class:: DihedralForceHarmonic(all_info)

   The constructor of dihedral harmonic interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string name, float k, float delta)
   
      specifies the dihedral harmonic force parameters with dihedral type, multiplicative constant, and phase shift angle.	

   .. py:function:: setCosFactor(float f)
   
      specifies the dihedral harmonic force parameters with factor.		  
	  
   Example::
   
      dihedralforce = gala.DihedralForceHarmonic(all_info)
      dihedralforce.setParams('A-B-B-A', 10.0, 0.0)
      app.add(dihedralforce)
	  
.. _harmonic-dihedral-improper:	
	  
Harmonic dihedral potential (improper) 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)=k\left( \varphi-\delta \right)^2		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - potential constant ``k`` (in units of energy/radians^2)
    - :math:`\delta` - phase shift angle ``delta`` (in radians)

    .. note::
	    The dihedral angles delta in script are in the unit of degree, and the program will convert them into radian automatically.	

.. py:class:: DihedralForceHarmonic(all_info)

   The constructor of dihedral harmonic interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string name, float k, float delta)
   
      specifies the dihedral harmonic force parameters with dihedral type, potential constant, and phase shift angle.	  
	  
   .. py:function:: setParams(string name, float k, float delta, int property)
   
      specifies the dihedral harmonic force parameters with dihedral type, potential constant, phase shift angle, and the property of proper or improper.	 	  
	  
   Example::
   
      dihedralforce = gala.DihedralForceHarmonic(all_info)
      dihedralforce.setParams('A-B-B-A', 10.0, 0.0, gala.DihedralForceHarmonic.Prop.improper)
      app.add(dihedralforce)

.. _opls-dihedral:	  
	  
OPLS dihedral potential (proper)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)=k_{1}+k_{2}\left[ 1+\cos \left( \varphi-\delta \right) \right]+k_{3} \left[ 1-\cos \left( 2\varphi-2\delta \right) \right] \\
		    + k_{4} \left[ 1+\cos \left( 3\varphi-3\delta \right) \right] + k_{5}\left[ 1 - \cos \left( 4\varphi-4\delta \right) \right]
        \end{eqnarray*}

    Coefficients:

    - :math:`k_1, k_2, k_3, k_4, k_5` - multiplicative constant ``k1, k2, k3, k4, k5`` (in units of energy)
    - :math:`\delta` - phase shift angle ``delta`` (in radians)
	
    .. note::
	    The dihedral angles delta in script are in the unit of degree, and the program will convert them into radian automatically.

.. py:class:: DihedralForceOPLSCosine(all_info)

   The constructor of dihedral OPLS cosine interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string name, float k1, float k2, float k3, float k4, float delta)
   
      specifies the dihedral OPLS cosine force parameters with dihedral type, k1, k2, k3, k4, and phase shift angle. In this function, the default k5 is 0.0.
	  
   .. py:function:: setParams(string name, float k1, float k2, float k3, float k4, loat k5, float delta)
   
      specifies the dihedral OPLS cosine force parameters with dihedral type, k1, k2, k3, k4, k5, and phase shift angle.	  
	  
   Example::
   
      dihedralforce = gala.DihedralForceOPLSCosine(all_info)
      dihedralforce.setParams('C_33-C_32-C_32-C_32', 0.0, 2.95188, -0.566963, 6.57940, 2.432826, 0.0)
      app.add(dihedralforce)
	  
.. _rb-dihedral:	
	  
Ryckaert-Bellemans potential (proper)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)= \sum_{n=0}^{5}C_{n}\left( \cos \left( \varphi \right) \right)^n
        \end{eqnarray*}

    Coefficients:

    - :math:`C_n` - multiplicative constant ``C0, C1, C2, C3, C4, C5`` (in units of energy)

.. py:class:: DihedralForceRyckaertBellemans(all_info)

   The constructor of dihedral RB interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string name, float c0, float c1, float c2, float c3, float c4, float c5)
   
      specifies the dihedral RB force parameters with dihedral type, c0, c1, c2, c3, c4, and c5.
	  
   Example::
   
      pfh = gala.DihedralForceRyckaertBellemans(all_info)
      pfh.setParams("A-A-B-B", 30.334, 0.0, -30.334, 0.0, 0.0, 0.0)
      app.add(pfh)

.. _amber-dihedral:	
	 
Amber potential (proper and improper)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{dihedral}}(\varphi)= k_{1} \left( 1 + \cos \left( \varphi - \delta_{1} \right) \right) +  k_{2} \left( 1 + \cos \left( 2\varphi - \delta_{2} \right) \right) \\
		+  k_{3} \left( 1 + \cos \left( 3\varphi - \delta_{3} \right) \right) +  k_{4} \left( 1 + \cos \left( 4\varphi - \delta_{4} \right) \right)
        \end{eqnarray*}

    Coefficients:

    - :math:`k_1, k_2, k_3, k_4` - multiplicative constant ``k1, k2, k3, k4`` (in units of energy)
    - :math:`\delta_1, \delta_2, \delta_3, \delta_4` - multiplicative constant ``delta1, delta2, delta3, delta4`` (in units of energy)	

    .. note::
	    The delta1, delta2, delta3, delta4 in script are in the unit of degree, and the program will convert them into radian automatically.

.. py:class:: DihedralForceAmberCosine(all_info)

   The constructor of dihedral Amber interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string name, float k1, float k2, float k3, float k4, float delta1, float delta2, float delta3, float delta4, int property)
   
      specifies the dihedral Amber force parameters with dihedral type, k1, k2, k3, k4, delta1, delta2, delta3, delta4, and the property of proper or improper.	
	  
   Example::
   
      pfh = gala.DihedralForceAmberCosine(all_info)
      pfh.setParams("A-A-B-B", 1.0, 0.0, 0.0, 0.0, 180.0, 0.0, 0.0, 0.0, 
		gala.DihedralForceAmberCosine.Prop.proper)
      app.add(pfh)	 


