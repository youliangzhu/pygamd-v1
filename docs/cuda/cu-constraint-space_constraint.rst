Space constraint
================

Bounce back condition
---------------------

.. py:class:: BounceBackConstrain(all_info, group)

   The constructor of a bounce back wall object with a group of particles.
   
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of charged particles.     
   

   .. py:function:: addWall(float o_x, float o_y, float o_z, float d_x, float d_y, float d_z)
   
      add wall with original point(o_x, o_y, o_z) and normal direction(d_x, d_y, d_z).
	  
   .. py:function:: addCylinder(float o_x, float o_y, float o_z, float d_x, float d_y, float d_z, float r)
   
      add cylinder with original point (o_x, o_y, o_z) , axis direction (d_x, d_y, d_z) and radius.
	  
   .. py:function:: addSphere(float o_x, float o_y, float o_z, float r)
   
      sphere with center point(o_x, o_y, o_z) and radius.
	  
   .. py:function:: clcearWall()
   
      clear the walls.
	  
   .. py:function:: clearCylinder()
   
      clear the cylinders
	  
   .. py:function:: clearSphere()
   
      clear the spheres.	  
	  
   
   Example::
   
      bbc = gala.BounceBackConstrain(all_info, group)
      bbc.addWall(0.0, 10.0, 0.0, 0.0, 1.0, 0.0)
      bbc.addWall(0.0, -10.0, 0.0, 0.0, 1.0, 0.0)
      app.add(bbc)


LJ surface force
----------------
 
.. py:class:: LJConstrainForce(all_info, group, r_cut)

   The constructor of a LJ interaction surface object for a group of particles.

   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of charged particles. 
   :param float r_cut: The cut-off radius.	   
   
   .. py:function:: setParams(string type, float epsilon, float sigma, float alpha)
   
      sets the interaction parameters of particle type for LJ surface.
	  
   .. py:function:: addWall(float o_x, float o_y, float o_z, float d_x, float d_y, float d_z)
   
      adds wall with original point(o_x, o_y, o_z) and normal direction(d_x, d_y, d_z) 
	  
   .. py:function:: addCylinder(float o_x, float o_y, float o_z, float d_x, float d_y, float d_z, float r)
   
      adds cylinder with original point(o_x, o_y, o_z) ,axis direction(d_x, d_y, d_z), and radius.
	  
   .. py:function:: addSphere(float o_x, float o_y, float o_z, float r)
   
      adds sphere with center point(o_x, o_y, o_z) and radius.
	  
   .. py:function:: clcearWall()
   
      clear the walls.
	  
   .. py:function:: clearCylinder()
   
      clear the cylinders
	  
   .. py:function:: clearSphere()
   
      clear the spheres.		  
	  
   Example::
   
      ljc = gala.LJConstrainForce(all_info, group, 1.0)
      ljc.addWall(0.0, 10.0, 0.0, 0.0, 1.0, 0.0)
      ljc.addWall(0.0, -10.0, 0.0, 0.0, 1.0, 0.0)
      ljc.setParams("A", 1.0, 1.0, 1.0)
      app.add(ljc)


