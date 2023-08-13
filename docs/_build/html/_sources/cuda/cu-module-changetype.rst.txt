Particle type change
====================

.. note::

A particle type change method can be used for the coarse-grained simulation of catalytic reaction of small molecules.

.. image:: changetype.png
    :width: 512 px
    :align: center
    :alt: illustration of particle type change method


ChangeType
----------

.. py:class:: ChangeType(all_info, source_type, target_type)

   The constructor of an object of ChangeType.
	 
   :param AllInfo all_info: The system information.
   :param string source_type: The type of source particles.  
   :param string target_type: The type of target particles. 

   .. py:function:: setPr(float prob)
   
      specifies reaction probability.
	  
   .. py:function:: setSite(NeighborList nlist, string site_type, float rcut)
   
      changing source particles to target particles, triggered by the sites (with site_type) in the neighbor list (nlist) of source particles within a cutoff of radius (rcut). 
	  
   .. py:function:: setWall(float o_x, float o_y, float o_z, float d_x, float d_y, float d_z)
   
      specifies the wall with original point(o_x, o_y, o_z) and normal direction(d_x, d_y, d_z), going through which the source particles will be changed to target particles.
	  
   .. py:function:: setInterface(NeighborList nlist, string type1, string type2, float rcut)
   
      specifies the interface between type1 particles and type2 particles, going through which the source particles will be changed to target particles. The interface is calculated from the neighbor list of sorce particles with a cutoff radius(rcut).

   .. py:function:: setNPTargetType(int np)
   
      specifies the number of target particles.
	  
   Example::
   
      ct = gala.ChangeType(all_info, "A", "B")
      ct.setSite(nlist, "C", 1.0)
      ct.setPr(0.002)
      ct.setPeriod(50)
      app.add(ct)
