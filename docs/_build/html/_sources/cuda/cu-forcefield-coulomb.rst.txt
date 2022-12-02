Coulomb interaction
===================

.. _ewald-theory:

Ewald summation theory
----------------------

The Coulomb interaction between two charge particles is given by:

   .. math::
       :nowrap:
   
       \begin{eqnarray*}
	   U\left( r \right)=f\frac{q_{i} q_{j}}{\epsilon_{r}r}
       \end{eqnarray*}
	   
where electric conversion factor :math:`f= 1/4\pi \epsilon_0=138.935\text{ }kJ\text{ }mol^{-1}\text{ }nm\text{ }e^{-2}`.
The total electrostatic energy of N particles and their periodic images is given by

   .. math::
       :nowrap:
   
       \begin{eqnarray*}
        V=\frac{f}{2\epsilon_{r}}\sum\limits_{\mathbf{n}}\sum\limits_{i}^{N}{\sum\limits_{j}^{N}{\frac{{q}_{i}{q}_{j}}{\left| {r}_{ij}+\mathbf{n} \right|}}}
       \end{eqnarray*}

The electrostatic potential is practically calculated by

   .. math::
       :nowrap:
   
       \begin{eqnarray*}
	   U\left( r^{*} \right)=\frac{q^{*}_{i} q^{*}_{j}}{r^{*}}
       \end{eqnarray*}
	   
The electric conversion factor and relative dielectric constant are considered in the reduced charge. 
For example, if the mass, length, and energy units are [amu], [nm], and [kJ/mol], respectively, according to :ref:`charge-units` the reduced charge is
:math:`q^{*}=z\sqrt{f^*/{\epsilon }_{r}}` with :math:`f^* = 138.935`. The :math:`z` is the valence of ion.

The calculation of Coulomb interaction is split into two parts, short-range part and long-range part by adding and subtracting a Gaussian distribution.

   .. math::
       :nowrap:
   
       \begin{eqnarray*}
	   G\left( r \right)=\frac{\kappa^{3}}{\pi^{3/2}}\mbox{exp}\left(-\kappa^2r^2\right)
       \end{eqnarray*}
	   
The short-range part including :py:class:`EwaldForce` and :py:class:`DPDEwaldForce` (for DPD) methods is calculated directly as non-bonded interactions.
The long-range part inlcuding :py:class:`PPPMForce` or :py:class:`ENUFForce` methods is calculated in the reciprocal sum by Fourier transform. 

For Coulomb interaction calculation, a short-range method and a long-range method are both needed.

   Example::
   
      groupC = gala.ParticleSet(all_info, "charge")
	  
      # real space
      ewald = gala.EwaldForce(all_info, neighbor_list, groupC, 3.0)#(,,r_cut)
      app.add(ewald)	  
	  
      # reciprocal space
      pppm = gala.PPPMForce(all_info, neighbor_list, groupC)
      pppm.setParams(32, 32, 32, 5, 3.0) 
      # grid number in x, y, and z directions, spread order, r_cut in real space.
      app.add(pppm)
      
      kappa = pppm.getKappa() 
      # an optimized kappa can be calculated by PPPMForce and passed into EwaldForce.
      ewald.setParams(kappa)

.. _ewald-short:	  
	  
Ewald (short-range)
-------------------------------------

Description:

    The short-range term is exactly handled in the direct sum.

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        V^{S}=\frac{f}{2\epsilon_{r}}\sum\limits_{\mathbf{n}}\sum\limits_{i}^{N}\sum\limits_{j}^{N}\frac{{q}_{i}{q}_{j}\mbox{erfc} \left(\kappa\left| {r}_{ij}+\mathbf{n} \right| \right)}{\left| {r}_{ij}+\mathbf{n} \right|}
       \end{eqnarray*}

    The following coefficients must be set:
	   
    - :math:`\kappa` - *kappa* (unitless)
	
.. py:class:: EwaldForce(all_info, nlist, group, r_cut)

   The constructor of an direct Ewald force object for a group of charged particles.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.
   :param ParticleSet group: The group of charged particles. 
   :param float r_cut: The cut-off radius.	  

   .. py:function:: setParams(string typei, string typej, float kappa)
   
      specifies the kappa per unique pair of particle types.
	  
   .. py:function:: setParams(float kappa)
   
      specifies the kappa for all pairs of particle types.
	  
   Example::
   
      group = gala.ParticleSet(all_info, "charge")
      kappa=0.8
      ewald = gala.EwaldForce(all_info, neighbor_list, group, 3.0)
      ewald.setParams(kappa)
      app.add(ewald)

.. _ewald-short-dpd:
	  
Ewald for DPD (short-range)
-------------------------------------

Description:

    In order to remove the divergency at :math:`r=0`, a Slater-type charge density is used to describe the charged DPD particles.

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        \rho(r)=\frac{q}{\pi\lambda^{3}}e^{-2r/\lambda}
       \end{eqnarray*}
	   
    - :math:`\lambda` - the decay length of the charge (in distance units)	
	
    The short-range term is exactly handled in the direct sum. 

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        V^{S}=\frac{f}{2\epsilon_{r}}\sum\limits_{\mathbf{n}}\sum\limits_{i}^{N}\sum\limits_{j}^{N}\frac{{q}_{i}{q}_{j}\mbox{erfc} \left(\kappa\left| {r}_{ij}+\mathbf{n} \right| \right)}{\left| {r}_{ij}+\mathbf{n} \right|} \left[1-(1+\beta r_{ij}\mbox{e}^{-2\beta r_{ij}} \right]
       \end{eqnarray*}

    The following coefficients must be set:
	   
    - :math:`\kappa` - *kappa* (unitless)
    - :math:`\beta=1/\lambda` - *beta* (in inverse distance units)	
	
.. py:class:: DPDEwaldForce(all_info, nlist, group, r_cut)

   The constructor of an direct Ewald force object for a group of charged particles.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.
   :param ParticleSet group: The group of charged particles. 
   :param float r_cut: The cut-off radius.	  

   .. py:function:: setParams(string typei, string typej, float kappa)
   
      specifies the kappa per unique pair of particle types.
	  
   .. py:function:: setParams(float kappa)
   
      specifies the kappa for all pairs of particle types.
	  
   .. py:function:: setBeta(float beta)
   
      specifies the beta for all pairs of particle types.  
	  
   Example::
   
      group = gala.ParticleSet(all_info, "charge")
      kappa=0.8
      dpd_ewald = gala.DPDEwaldForce(all_info, neighbor_list, group, 3.0)
      dpd_ewald.setParams(kappa)
      app.add(dpd_ewald)	  

.. _pppm-long:	  
	  
PPPM (long-range)
----------------------

Description:

    The long-range term is exactly handled in the reciprocal sum. 

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        V^{L}&=&\frac{1}{2V\epsilon_{0}\epsilon_{r}}\sum\limits_{\mathbf{k}\neq0}\frac{\mbox{exp}(-\mathbf{k}^{2}/4\kappa^{2})}{\mathbf{k}^{2}} \left| S(\mathbf{k}) \right|^{2} \\
        S(\mathbf{k})&=&\sum\limits_{i=1}^{N}q_{i}\mbox{exp}^{i\mathbf{k} \cdot \mathbf{r}_i}		
       \end{eqnarray*}

    The self-energy term. 

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        V^{self}&=&\frac{1}{f}\frac{\kappa}{\sqrt{\pi}}\sum\limits_{i=1}^{N}q_{i}^{2}		
       \end{eqnarray*}	
	   
    - :math:`\kappa` - *kappa* (unitless)
	
.. py:class:: PPPMForce(all_info, nlist, group)
	  
   The constructor of a PPPM force object for a group of charged particles.

   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.
   :param ParticleSet group: The group of charged particles.

   .. py:function:: setParams(int nx, int ny, int nz, int order, float r_cut)
   
      specifies the PPPM force with the number of grid points in x, y, and z direction, the order of interpolation, and the cutoff radius of direct force.
	  
   .. py:function:: setParams(float fourierspace, int order, float r_cut)
   
      specifies the PPPM force with the fourier space, the order of interpolation, and the cutoff radius of direct force.
      The number of grid points will be derived automatically.
	  
   .. py:function:: float getKappa()
   
      return the kappa calculated by PPPM force.
	  
   Example::
   
      group = gala.ParticleSet(all_info, "charge")
      pppm = gala.PPPMForce(all_info, neighbor_list, group)
      pppm.setParams(32, 32, 32, 5, 3.0)
      app.add(pppm)

.. _enuf-long:
	  
ENUF (long-range)
----------------------

.. py:class:: ENUFForce(all_info, nlist, group)
	  
   The constructor of an ENUF force object for a group of charged particles.

   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.
   :param ParticleSet group: The group of charged particles.

   .. py:function:: setParams(float alpha, float sigma, int precision, int Nx, int Ny, int Nz)
      
      specifies the ENUF force with alpha, hyper sampling factor sigma, precision determine the order of interpolation (precision*2+2), and the number of grid points in x, y, and z direction.	
	  
   Example::
   
      group = gala.ParticleSet(all_info, "charge")
      kappa=0.8
      enuf = gala.ENUFForce(all_info, neighbor_list, group)
      enuf.setParams(kappa, 2.0, 2, 32, 32, 32)
      app.add(enuf)

	   