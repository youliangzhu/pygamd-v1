External interaction
====================

Axial stretching
----------------

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


External force
--------------

.. py:class:: ExternalForce(all_info, group)

   The constructor of an external force object for a group of particles. The external force will be added on 
   each particle.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.   

   .. py:function:: setForce(Variant vf, std::string direction)
   
      specifies the force magnitude varying by time steps and direction (the candidates are "X", "Y", and "Z").
   
   .. py:function:: setForce(Variant vf, float x, float y, float z)
   
      specifies the force magnitude varying by time steps and direction vector (x, y, z).
   
   .. py:function:: setParams(string type, float factor)
   
      specifies the factor of external force for a particle type (the default value is 1.0).
   
   .. py:function:: setParams(unsigned int index, float factor)
   
      specifies the factor of external force for a particle with index.
	  
   Example::
   
      v = gala.VariantSin()
      v.setPoint(0, 1000, 1, -1) 
      v.setPoint(1000000, 1000, 1, -1)
      # set the parameters of sinusoid force by time step, period, max and min value where 
	  # the latter three parameters are linearly varying by time step.
	  
      groupA = gala.ParticleSet(all_info, "A")
      ef = gala.ExternalForce(all_info, groupA)
      #initializes an external force object with system information and particle group.
      ef.setForce(v, "X")
      # sets parameters with force and direction.
      app.add(ef)



