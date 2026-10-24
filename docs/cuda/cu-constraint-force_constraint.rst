Force constraint
================

Fixed force harmonic
--------------------
 
.. py:class:: FixedForceHarmonic(all_info, group)

   The constructor of a fixed harmonic force object for a group of particles. 
   The fixed harmonic force can be applied to the particles at X, Y, Z, XY, YZ, XZ, and XYZ by the free combination of setFixedPosition() functions.

   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of constraint particles.
   
   .. py:function:: setParams(string type, float K, float r_0)
   
      sets the harmonic interaction parameters for the type particles where K is the strength and r_0 is the equilibrium length.
	  
   .. py:function:: setFixedPositionX(float x)
   
      set the fixed position at X direction. 
	  
   .. py:function:: setFixedPositionY(float y)
   
      set the fixed position at Y direction. 
	  
   .. py:function:: setFixedPositionZ(float z)
   
      set the fixed position at Z direction. 

   Example::
   
      ffh = gala.FixedForceHarmonic(all_info, group)
      ffh.setParams("A", 100.0, 1.0) #particle type, harmonic constant, equilibrium length
      ffh.setFixedPositionZ(0.0) #fixed position at Z direction
      app.add(ffh)


