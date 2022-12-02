Dissipative particle dynamics
=============================

DPD force
---------

Description:

    The DPD force consists of pair‐wise conservative, dissipative and random terms.

    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        \vec{F}_{ij}^{C}&=&\alpha\left(1-\frac{r_{ij}}{r_{cut}}\right)\vec{e}_{ij} \\
        \vec{F}_{ij}^{D}&=&-\gamma\omega^{D}(r_{ij})(\vec{e}_{ij} \cdot \vec{v}_{ij} )\vec{e}_{ij}  \\	
        \vec{F}_{ij}^{R}&=&T\sigma\omega^{R}(r_{ij})\xi_{ij}\vec{e}_{ij} \\			
       \end{eqnarray*}

	   
    - :math:`\gamma=\sigma^{2}/2k_{B}T`
    - :math:`\omega^{D}(r_{ij})=[\omega^{R}(r_{ij})]^2=(1-r_{ij}/r_{cut})^2`	
    - :math:`\xi_{ij}` - a random number with zero mean and unit variance
    - :math:`T` - `temperature`
      - *optional*: defaults to 1.0	
    - :math:`r_{cut}` - *r_cut* (in distance units)	
      - *optional*: defaults to 1.0

    The following coefficients must be set per unique pair of particle types:
	
    - :math:`\alpha` - *alpha* (in energy units)
    - :math:`\sigma` - *sigma* (unitless)


.. py:class:: DPDForce(all_info, nlist, r_cut, temperature, rand_num)

   The constructor of a DPD interaction object.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.
   :param float temperature: The temperature.
   :param int rand_num: The seed of random number generator.   
   
.. py:class:: DPDForce(all_info, nlist, r_cut, rand_num)

   The constructor of a DPD interaction object. The default temperature is 1.0.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.
   :param int rand_num: The seed of random number generator.
   
   .. py:function:: setParams(string typei, string typej, float alpha, float sigma)
   
      specifies the DPD interaction parameters per unique pair of particle types.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time step.
	  
   .. py:function:: setDPDVV()
   
      calls the function to enable DPDVV method (the default is GWVV).
	  
   Example::
   
      dpd = gala.DPDForce(all_info, neighbor_list, 1.0, 12345)
      dpd.setParams('A', 'A', 25.0, 3.0)
      app.add(dpd)
	  
GWVV integration
----------------

Description:

    Integration algorithm.


    .. math::
       :nowrap:
   
       \begin{eqnarray*}
        &v_i^0&\leftarrow v_i+ \lambda\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})\\ 
        &v_i&\leftarrow v_i+ \frac{1}{2}\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})\\
        &r_i&\leftarrow r_i+ v_i \Delta t\\
        &&Calculate \quad F_i^c\{r_j\}, F_i^d\{r_j, v_j^0\}, F_i^r\{r_j\}\\
        &v_i&\leftarrow v_i+ \frac{1}{2}\frac{1}{m}(F_i^c \Delta t + F_i^d \Delta t + F_i^r \sqrt{\Delta t})		
       \end{eqnarray*}

    - :math:`\lambda` - *lambda* (unitless)
      - *optional*: defaults to 0.65
	  
.. py:class:: DPDGWVV(AllInfo all_info, ParticleSet group)

   The constructor of a GWVV NVT thermostat for a group of DPD particles.

   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	   

   .. py:function:: setLambda(float lambda)
   
      specifies lambda.	  
	  
   Example::

      gwvv = gala.DPDGWVV(all_info, group)
      app.add(gwvv)
	  
Coulomb interaction in DPD
--------------------------

Description:

    In order to remove the divergency at :math:`r=0`, a Slater-type charge density is used to describe the charged DPD particles.
    Thereby, :ref:`ewald-short-dpd` (:py:class:`DPDEwaldForce`) method can be employed to calculated the short-range part of Ewald summation.
    The long-range part of Ewald summation can be calculated by :ref:`pppm-long` (:py:class:`PPPMForce`) or :ref:`enuf-long` (:py:class:`ENUFForce`).
    And the :ref:`enuf-long` (:py:class:`ENUFForce`) is suggested.
 
   Example::

      group_charge = gala.ParticleSet(all_info, "charge")
      kappa=0.2

      # real space
      ewald = gala.DPDEwaldForce(all_info, neighbor_list, group_charge, 3.64)#(,,rcut)
      ewald.setParams(kappa)
      app.add(ewald)
	  
      # reciprocal space 
      enuf = gala.ENUFForce(all_info, neighbor_list, group_charge)
      enuf.setParams(kappa, 2, 2, 20, 20, 20)
      app.add(enuf)
 
 
Reduced charges:

    The charges should be converted into the ones in reduced units according to :ref:`charge-units`.
    Typically, the fundamental length and energy are :math:`\sigma=0.646\text{ }nm` and :math:`\epsilon=k_{B}T` with :math:`T=300\text{ }K`, respectively, in DPD.
    The reduced charges are :math:`q^{*}=z\sqrt{f/(\sigma k_{B}T \epsilon_r)}`. The :math:`z` is the valence of ion.

    Here is a :ref:`molgen` script for polyelectrolyte. 

   Example::
   
      #!/usr/bin/python
      import sys

      import molgen
      import math
      
      er=78.0
      kBT=300.0*8.314/1000.0
      r=0.646
      gama=138.935
      dpdcharge=math.sqrt(gama/(er*kBT*r))
      
      mol1=molgen.Molecule(50)
      mol1.setParticleTypes("P*50")
      topo="0-1"
      for i in range(1,50-1):
      	c=","+str(i)+"-"+str(i+1)
      	topo+=c
      mol1.setTopology(topo)
      mol1.setBondLength(0.7)
      mol1.setMass(1.0)
      
      mol2=molgen.Molecule(1)
      mol2.setParticleTypes("C")
      mol2.setMass(1.0)
      mol2.setCharge(dpdcharge)
      
      mol3=molgen.Molecule(1)
      mol3.setParticleTypes("A")
      mol3.setMass(1.0)
      mol3.setCharge(-dpdcharge)
      
      mol4=molgen.Molecule(1)
      mol4.setParticleTypes("W")
      mol4.setMass(1.0)
      
      gen=molgen.Generators(15,15,15)
      gen.addMolecule(mol1,1)
      gen.addMolecule(mol2,75)
      gen.addMolecule(mol3,75)
      gen.addMolecule(mol4,9925)
      
      gen.outPutXml("ps0")	  
