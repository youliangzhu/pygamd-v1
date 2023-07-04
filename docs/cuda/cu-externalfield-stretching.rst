Axial stretching
================


.. py:class:: AxialStretching(all_info, group)

   The constructor of a stretching object of the box for a group of particles.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.
   

   .. py:function:: setBoxLength(Variant vL, string direction)
   
      specifies the change of box length and its direction with Variant. The candidates are "X", "Y", "Z".
	  
   Example::
   
      v = gala.VariantLinear()
      v.setPoint(0, 31) # time step, box length.
      v.setPoint(100000, 60)
	  
      axs = gala.AxialStretching(all_info, group)
      axs.setBoxLength(v, 'Y')
      app.add(axs)


