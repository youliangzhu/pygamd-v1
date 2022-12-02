Anisotropic particle
====================

Gay-Berne model
---------------

Uniaxial GB interaction
^^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: GBForce(all_info, nlist, r_cut)

   The constructor of a method of Gay-Berne force.
   
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.	  

   .. py:function:: setParams(string type1, string type2, float epsilon0, float sigma0, float nu, float mu, float sigma_e, float sigma_s, float epsilon_e, float epsilon_s, float Ps)
   
      specifies the GB force parameters with type1, type2, epsilon0, sigma0, nu, mu, end-to-end length (sigma_e), side-by-side length (sigma_s), end-to-end energy (epsilon_e), side-by-side energy (epsilon_s), Ps.
	  
   Example::
   
      gb = gala.GBForce(all_info, neighbor_list, 10.0)
      gb.setParams('A', 'A' , 1.5, 1.5, 1.0, 2.0,3.0, 1.0, 0.5, 3.0, 1.0)
      # sets parameters: type1, type2, epsilon0, sigma0, nu, mu, sigma_e, 
      # sigma_s, epsilon_e, epsilon_s, Ps.
      app.add(gb)
	  
Bond force of uniaxial GB particles
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: BondForceAni(all_info)

   The constructor of a method of bond force calculation of anisotropic particles.
   
   :param AllInfo all_info: The system information.   

   .. py:function:: setParams(string bondtype, float Kbond, float rbond, float Kangle, float dangle)
   
      specifies the bond force parameters with bond type, bond spring constant, end-to-end length of GB particle, angle spring constant, equilibrium angle degree.
	  
   Example::
   
      bondani = gala.BondForceAni(all_info)
      bondani.setParams('A-A', 100.0 , 4.498, 30.0, 0.0)
      app.add(bondani)
	  
	  
Biaxial GB interaction
^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: PBGBForce(all_info, nlist)

   The constructor of a method of Gay-Berne force.
   
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.   

   .. py:function:: setGUM(float gamma, float nu, float mu)
   
      specifies the GB force parameters with gamma, nu, mu.

   .. py:function:: setParams(string type1, string type2, float epsilon, float sigma, float r_cut)
   
      specifies the GB force parameters with type1, type2, epsilon, sigma, cutoff radius.
	  
   .. py:function:: setAspheres(string filename)
   
      specifies the file for shape parameters.

   .. py:function:: setPatches(string filename)
   
      specifies the file for Patch parameters.		  
	 
	  
   Example::
   
		pbgb = gala.PBGBForce(all_info, neighbor_list)
		pbgb.setGUM(1.0, 3.0, 1.0);#(gamma, niu, miu)
		pbgb.setParams('B', 'B' , 1.0, 1.0, 5.0)#(,,epsilon, sigma, rcut)
		pbgb.setParams('A', 'A' , 1.0, 1.0, 5.0)#(,,epsilon, sigma, rcut)
		pbgb.setParams('A', 'B' , 2.0, 1.0, 5.0)#(,,epsilon, sigma, rcut)
		pbgb.setAspheres('patch.log')#(,a,b,c,eia_one,eib_one,eic_one)
		pbgb.setPatches('patch.log')
		app.add(pbgb)

   File 'patch.log'::
	
		<Patches>
		B 2                             #particle type, patch number
		p1 60  0         0   1          #patch type, beta(degree) which is half of the opening angle-
		p1 60  0         0   -1         #of the attractive patch, patch position(x, y, z) in unit vector
		</Patches>
		<PatchParams>
		p1 p1 88.0 0.5                  #patch type, patch type, alpha_A, and gamma_epsilon
		</PatchParams>
		<Aspheres>
		A 1.0 1.0 1.0 3.0 3.0 3.0       #a,b,c,eia_one,eib_one,eic_one
		B 1.0 1.0 3.0 1.0 1.0 0.2       #a,b,c,eia_one,eib_one,eic_one
		</Aspheres>  

Soft anisotropic model
----------------------

Janus particle model
^^^^^^^^^^^^^^^^^^^^

.. py:class:: LZWForce(all_info, nlist, r_cut)

   The constructor of a method of LZW force calculation.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.	  

   .. py:function:: setParams(string type1, string type2, float alphaR, float mu, float nu, float alphaA, float beta)
   
      specifies the LZW force parameters with type1, type2, alphaR, mu, nu, alphaA, and beta.
	  
   .. py:function:: setMethod(string method)
   
      chooses a method of 'Disk', 'Janus', ABAtriJanus', 'BABtriJanus'.
	  
   Example::
   
      lzw = gala.LZWForce(all_info, neighbor_list, 1.0)
      lzw.setParams('A', 'A' , 396.0, 1.0, 0.5, 88.0,60.0/180.0*3.1415926)
      lzw.setMethod('ABAtriJanus')
      # sets method with the choice of ABAtriJanus.
      app.add(lzw)

Thermostat for Janus particle model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. py:class:: BerendsenAniNVT(AllInfo all_info, ParticleSet group, ComputeInfo comp_info, float T, float tauT, float tauR)

   The constructor of a Berendsen NVT thermostat for anisotropic particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param ComputeInfo comp_info: The object of calculation of collective information.	   
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.	 	  

   .. py:function:: setTau(float tauT, float tauR)
   
      specifies the Berendsen NVT thermostat with tauT and tauR.
	  
   .. py:function:: setT(float T)
   
      specifies the temperature with a constant value.
	  
   .. py:function:: setT(Variant vT)
   
      specifies the temperature with a varying value by time step.
	  
   Example::
   
      bere = gala.BerendsenAniNVT(all_info, group, comp_info, 1.0, 0.3, 0.1)
      app.add(bere)
	  
Multiple patch particle model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Description:

    .. math::
        :nowrap:

        \begin{eqnarray*}
		U_{ij}=\left\{ \begin{array}{ll} \frac{\alpha_{ij}^R d_{ij}}{\mu}\left(1-\frac{r_{ij}}{d_{ij}}\right)^\mu- 
		\sum\limits_{\kappa=1}^{M_i}\sum\limits_{\lambda=1}^{M_j} f^\nu \left(\mathbf{n}_{i}^{\kappa}, \mathbf{n}_j^{\lambda},  \mathbf{r}_{ij}\right) \frac{\alpha_{ij}^A d_{ij}}{\mu}\left[\frac{r_{ij}}{d_{ij}}-\left(\frac{r_{ij}}{d_{ij}}\right)^\mu\right] &  r_{ij}\leq d_{ij} \\  0  & r_{ij}> d_{ij}, 
		\end{array} \right. 
        \end{eqnarray*}
		
		\begin{eqnarray*} 
		f\left(\mathbf{n}_{i}^{\kappa}, \mathbf{n}_j^{\lambda},  \mathbf{r}_{ij}\right) = \left\{ \begin{array}{ll} 
		\cos\frac{\pi\theta_i^{\kappa}}{2\theta_{m}^{\kappa}}\cos\frac{\pi\theta_j^{\lambda}}{2\theta_{m}^{\lambda}} & \textrm{if $\cos\theta_i^{\kappa}\geq \cos\theta_{m}^{\kappa}$ and $\cos\theta_j^{\lambda} \geq \cos\theta_{m}^{\lambda}$}\\ 
		0 & \textrm{otherwise}. 
		\end{array} \right. 
		\end{eqnarray*}		

    The following coefficients must be set per unique pair of particle types:

    - :math:`\alpha^R` - *alphaR*, repulsive strength
    - :math:`\mu` - *mu*, the power (unitless)
    - :math:`\alpha^A` - *alphaA*, attractive strength
    - :math:`d` - the diameter defaults to the r_cut (in distance units)
    - :math:`\nu` - *nu*, the angular width of attraction (unitless)

.. py:class:: AniForce(all_info, nlist, r_cut)

   The constructor of force calculation of multiple patch particle model.
	  
   :param AllInfo all_info: The system information.
   :param NeighborList nlist: The neighbor list.  
   :param float r_cut: The cut-off radius.	  

   .. py:function:: setParams(string type1, string type2, float alphaR, float mu)
   
      specifies the force parameters with type1, type2, alphaR, mu.
	  
   .. py:function:: setPatches(string filename)
   
      specifies the file for Patch parameters.
	  
   Example::
   
		ani = gala.AniForce(all_info, neighbor_list, 1.0)
		ani.setParams('A', 'A' , 396.0, 2.0)#(,,alpha_R,mu)
		ani.setParams('A', 'B' , 396.0, 2.0)#(,,alpha_R,mu)
		ani.setParams('B', 'B' , 396.0, 2.0)#(,,alpha_R,mu)
		ani.setPatches('patch-3.log')
		app.add(ani)
		
   File 'patch-3.log'::
	
		<Patches>                       
		A 0                             #particle type, patch number
		B 3                             #particle type, patch number
		p1 45  0         0   1          #patch type, beta(degree) which is half of the opening angle-
		p2 45  0.866025  0   -0.5       #of the attractive patch, patch position(x, y, z) in unit vector
		p3 45 -0.866025  0   -0.5       
		</Patches>
		<PatchParams>
		p1 p1 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		p2 p2 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		p3 p3 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		p1 p2 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		p1 p3 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		p2 p3 220.0 0.5                 #patch type, patch type, alpha_A, and nu
		</PatchParams>	

Thermostat for multiple patch particle model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

   The motion of anisotropic particles with multiple patches are integrated by rigid body method with a body index in XML file. 
   The solvent particles with a body index of -1 are integrated by normal methods.

	  
   Example::
   
		bgroup = gala.ParticleSet(all_info, 'body')
		rigidnvt = gala.NVTRigid(all_info, bgroup, 1.0, 0.2)
		app.add(rigidnvt)
		
		nbgroup = gala.ParticleSet(all_info,'non_body')
		comp_info_nb = gala.ComputeInfo(all_info, nbgroup)
		nh = gala.NoseHooverNVT(all_info, nbgroup, comp_info_nb, 1.0, 1.0)#( ,temperature, tau)
		app.add(nh)
