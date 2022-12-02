Bond stretching
---------------

**Overview**

Bonds impose connected forces on specific pairs of particles to model chemical bonds.
The bonds are specified in :ref:`xml-format` configuration file with the format::

   <bond>
   bond_type(str)  particle_i(int)  particle_j(int)
   ...
   </bond>
   
By themselves, bonds do nothing. Only when you specify a bond force in script(i.e. :py:class:`BondForceHarmonic`), are forces actually calculated between the listed particles.

=======================   ===============================
:ref:`harmonic-bond`      :py:class:`BondForceHarmonic`
:ref:`fene-bond`          :py:class:`BondForceFene`
:ref:`polynominal-bond`   :py:class:`BondForcePolynomial`
:ref:`morse-bond`         :py:class:`BondForceMorse`
=======================   ===============================

.. image:: bond.png
    :width: 250 px
    :align: center
    :alt: Principle of bond stretching

.. _harmonic-bond:
	
Harmonic bond potential
^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{bond}}(r) = \frac{1}{2}k\left( r - r_{0} \right)^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - spring constant ``k`` (in units of energy/distance^2)
    - :math:`r_0` - equilibrium length ``r0`` (in distance units)

.. py:class:: BondForceHarmonic(all_info)

   The constructor of harmonic bond interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string type, float k, float r0)
   
      specifies the bond interaction parameters with bond type, spring constant, and equilibrium length.

   Example::
   
      bondforce = gala.BondForceHarmonic(all_info)
      bondforce.setParams('polymer', 1250.000, 0.470)
      app.add(bondforce)

.. _fene-bond:	  
	  
FENE bond potential
^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{bond}}(r)=-\frac{1}{2}k {r_{m}}^{2}\log \left[ 1-\frac{\left(r - r_{0} -\Delta \right)^{2}}{r_{m}^{2}} \right]
        \end{eqnarray*}

        \begin{eqnarray*}
           V_{\mathrm{WCA}}(r)=&4 \epsilon \left[ \left( \frac{\sigma }{r-\Delta } \right)^{12}-\left( \frac{\sigma }{r-\Delta } \right)^{6} \right] + \epsilon 
		                       & ,(r - \Delta)<\sigma^{1/6}  \\
                            = & 0 & ,(r - \Delta) \ge \sigma^{1/6}  \\
        \end{eqnarray*}
		
    Coefficients:

    - :math:`k` - attractive force strength ``k`` (in units of energy/distance^2)
    - :math:`r_0` - equilibrium length ``r0`` (in distance units)
      - *optional*: defaults to 0.0
    - :math:`r_m` - maximum bond length ``rm`` (in distance units)
    - :math:`\epsilon` - *epsilon* (in energy units)
    - :math:`\sigma` - *sigma* (in distance units)
    - :math:`\Delta = (d_{i} + d_{j})/2 - 1.0` - (in distance units); :math:`d_{i}` and :math:`d_{j}` are the diameter of particle :math:`i` and :math:`j` which can be input from XML file.	
    
    Note:
        :math:`\Delta` only will be considered (default value is 0.0) by calling the function ``setConsiderDiameter(True)`` 	
	
.. py:class:: BondForceFENE(all_info)

   The constructor of FENE bond interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string type, float k, float rm)
   
      specifies the FENE bond force parameters with bond type, spring constant, and the maximum length of the bond.
	  
   .. py:function:: setParams(string type, float k, float rm, float r0)
   
      specifies the FENE bond force parameters with bond type, spring constant, maximum length, and equilibrium length.
	  
   .. py:function:: setParams(string type, float k, float rm, float epsilon, float sigma)
   
      specifies the FENE+WCA bond parameters with bond type, spring constant, maximum length of the bond, epsilon, sigma (the latter two parameters for WCA force between two bonded particles ).

   .. py:function:: setConsiderDiameter(bool con_dia)
   
      the diameter of particles will be considered or not

   Example::
   
      bondforcefene = gala.BondForceFENE(all_info)
      bondforcefene.setParams('polymer', 10, 1.2)
      app.add(bondforcefene)

.. _polynominal-bond:		  
	  
Polynominal bond potential
^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{bond}}(r)=k_{1}\left( r - r_{0} \right)^{2}+k_{2}\left( r - r_{0} \right)^{4}
        \end{eqnarray*}

    Coefficients:

    - :math:`k_1` - spring constant ``k1`` (in units of energy/distance^2)
    - :math:`k_2` - spring constant ``k2`` (in units of energy/distance^4)	
    - :math:`r_0` - equilibrium length ``r0`` (in distance units)

	
.. py:class:: BondForcePolynomial(all_info)

   The constructor of polynomial bond interaction object.
 
   :param AllInfo all_info: The system information.

   .. py:function:: setParams(string type, float k1, float k2, float r0)
   
      specifies the polynomial bond force parameters with bond type, spring constant k1, spring constant k2, and equilibrium bond length r0.
	  
   Example::
   
      bondforce_polynomial = gala.BondForcePolynomial(all_info)
      bondforce_polynomial.setParams('polymer', 10.0, 100.0, 1.2)
      app.add(bondforce_polynomial)

.. _morse-bond:	  
	
Morse bond potential
^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        V_{\mathrm{bond}}(r)=&k\left[ 1-e^{-\alpha \left( r-r_{0} \right)} \right]^{2} & r < r_{\mathrm{m}} \\
                            = & 0 & r \ge r_{\mathrm{m}} \\		
        \end{eqnarray*}

    Coefficients:

    - :math:`k` - well depth ``k`` (in units of energy)
    - :math:`\alpha` - controls the 'width' of the potential ``alpha`` (he smaller :math:`\alpha` is, the larger the well)	
    - :math:`r_0` - equilibrium length ``r0`` (in distance units)
    - :math:`r_m` - maximum interaction range ``rm`` (in distance units)
	
.. py:class:: BondForceMorse(all_info)

   The constructor of Morse bond interaction object.
 
   :param AllInfo all_info: The system information.
   
   .. py:function:: setParams(string name, float k, float alpha, float r0, float rm)
   
      specifies the Morse bond force parameters with bond type, spring constant, alpha controls the 'width' of the potential, equilibrium bond length, maximum interaction range.
	  
   Example::
   
      bondforce_morse = gala.BondForceMorse(all_info)
      bondforce_morse.setParams('polymer', 10.0, 1.0, 1.0, 2.0)
      app.add(bondforce_morse)

	  
