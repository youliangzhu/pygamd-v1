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
        \vec{F}_{ij}^{R}&=&\sigma\omega^{R}(r_{ij})\xi_{ij}\vec{e}_{ij} \\			
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


.. py:class:: force.dpd(info, rcut=1.0)

   Constructor of a DPD interaction object.
	  
   :param info: system information.
   :param rcut: the cut-off radius of interactions. 

   .. py:function:: setParams(type_i, type_j, alpha, sigma)
   
      specifies the DPD interaction parameters between two types of particles.
	  
   Example::
   
      fn = pygamd.force.dpd(info=mst, rcut=1.0)
      fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
      app.add(fn)
	  
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
	  
.. py:class:: integration.gwvv(info, group)

   Constructor of a GWVV NVT thermostat for a group of DPD particles.

   :param info: system information.
   :param group: a group of particles.
	  
   Example::

      gw = pygamd.integration.gwvv(info=mst, group='all')
      app.add(gw)
	  
