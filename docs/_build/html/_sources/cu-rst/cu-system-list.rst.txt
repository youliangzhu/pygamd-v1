Particle list
=============

Neighbor list
-------------

.. py:class:: NeighborList(all_info, r_cut, r_buffer)

   Constructor of a neighbor list object
  
   :param AllInfo all_info: System information
   :param float r_cut: Cut-off radius
   :param float r_buffer: Buffer distance

   .. py:function:: setRCut(float r_cut, float r_buffer)
   
      specifies the cut-off and buffer distance.

   .. py:function:: setRCutPair(string typi, string typj, float r_cut) 

      specifies the cut-off per unique pair of particle types.
	  
   .. py:function:: setNsq()
   
      switches on the method of searching all particle to build up list.
	  
   .. py:function:: setDataReproducibility()
   
      switches on the data reproducibility.
	  
   .. py:function:: addExclusionsFromBonds()
   
      adds 1-2 exclusion into exclusion list.
	  
   .. py:function:: addExclusionsFromAngles()
   
      adds 1-3 exclusion into exclusion list.
	  
   .. py:function:: addExclusionsFromDihedrals()
   
      adds 1-4 exclusion into exclusion list.
	  
   .. py:function:: addExclusionsFromBodys()
   
      adds body exclusion into exclusion list.
	  
   .. py:function:: setFilterDiameters()
   
      considers the radius of particle in neighbor list which includes the particles within r_cut + (diameter_i + diameter_j)/2.
   
   Example::
   
      neighbor_list = gala.NeighborList(all_info, 3.0 ,0.4)

Cell list
---------

.. py:class:: CellList(all_info)

   Constructor of a cell list object

   :param AllInfo all_info: System information

   .. py:function:: setNominalWidth(float width)
   
      specifies the length of cell.
	  
   .. py:function:: setNominalDim(unsigned int x, unsigned int y, unsigned int z)
   
      specifies the dimensions of grid in 'X', 'Y', and 'Z' directions.
	  
   .. py:function:: setDataReproducibility()
   
      switches on data reproducibility function.
	  
   Example::
   
      cell_list = gala.CellList(all_info)

