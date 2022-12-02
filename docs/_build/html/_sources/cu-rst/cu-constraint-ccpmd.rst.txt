Constant chemical potential
===========================

.. note::

Constant Chemical Potential Molecular Dynamics (C :math:`\mu` MD) method introduces an external force that controls the environment of the chemical process of interest. 
This external force, drawing molecules from a finite reservoir, maintains the chemical potential constant in the region where the process takes place. 
This method is able to study crystal growth dynamics under constant supersaturation conditions or evaporation dynamics under constant vapor pressure. Reference: C. Perego, M. Salvalaglio, and M. Parrinello, J. Chem. Phys., 2015, 142, 144113.

.. image:: ccpmd.png
    :width: 512 px
    :align: center
    :alt: Constant Chemical Potential Molecular Dynamics

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
        F^{\mathrm{\mu}}(z) &=& k(n^{CR}-n_{0})G(z, Z_{F}) \\
		G(z-Z_{F}) &=& \frac{1}{4\omega}\left[ 1+cosh \left(\frac{z-Z_{F}}{\omega} \right) \right]^{-1} \\
        \end{eqnarray*}

    Coefficients:

    - :math:`n_0` - target constant concentration ``n0`` (in reduced units)
    - :math:`k` - spring constant ``k`` (in units of energy/distance^2)
    - :math:`\sigma` - width of external force region ``sigma`` (in units of distance)	
    - :math:`\omega` - an intensity peak proportional to ``1/omega`` and a width proportional to ``omega`` (in units of distance)	

.. py:class:: CCPMD(all_info, group)

   The constructor of a constant chemical potential object of a group of particles.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.
  
   .. py:function:: setWall(float ox, float oy, float oz, float dx, float dy, float dz)
   
      specifies external force region with plane center (ox, oy, oz) and direction (dx, dy, dz). If the normal direction of wall is in Z direction, the center position of plane is (0.0, 0.0, :math:`Z_{F}`).

   .. py:function:: setParams(float n0, float k, float sigma, float omega, float Len_CR)
   
      specifies target concentration, spring constant, sigma, omega, and the length of control region.
	  
   Example::

      groupS = gala.ParticleSet(all_info, 'S')
      ccp = gala.CCPMD(all_info, groupS)
      ccp.setWall(0.0, 0.0, -25.0, 0.0, 0.0, -1.0)
      ccp.setParams(0.5, 1000.0, 1.0, 0.1, 5.0)
      app.add(ccp)
  



