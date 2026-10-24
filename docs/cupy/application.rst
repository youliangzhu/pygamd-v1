Application
===========

Modules management
------------------

PYGAMD is organized by being composed of modules. Application manages and calls modules, and thereby run simulations.
Usually, we only define an application object in the context of script. The modules can be 
added into by ``add()`` or removed from by ``remove()`` the application before runing ``run()`` the simulation. 

.. py:class:: application.dynamics(info, dt, sort=True)

   Constructor of application object.
	  
   :param info: system information
   :param dt: integration time step
   :param sort: if device memory is sorted by Hilbert curve, the default is True. 	    

   .. py:function:: add(object)
   
      adds an object to the application.
	  
   .. py:function:: remove(object)
   
      removes an added object.
	
   .. py:function:: run(N)
   
      runs the simulation for N time steps.
	  
   Example::
 
      app = pygamd.application.dynamics(info=mst, dt=0.001)
      # builds up an application.
      app.run(10000)
      # runs the simulation for 10000 time steps.


Multi-stage simulation 
----------------------

An application can have single or multiple stage simulations. The commands in the context of script are executed sequentially.
Every stage simulation is achieved with ``run()``.  Before a stage of simulation, the modules and parameters can be adjusted.
New modules can be added into the applicaitons by ``add()``. The added modules at last stage can be removed from the application, 
otherwise they will be kept. For example:

* First stage simulation::

   app = pygamd.application.dynamics(info=mst, dt=0.001)
   app.add(lj)
   app.add(nvt)
   app.run(1000)
  
* Second state simulation::

   app.remove(lj)
   app.remove(nvt)
   app.add(harmonic)
   app.add(npt)
   app.run(1000)
   
   
Two-dimensional simulation 
--------------------------

1. Controlling script, i.e. 'file.py' script is same for two- and three-dimensional simulations.
 
2. However, configuration file i.e. MST file should indicate two-dimensional 
system by:

   1. pointing out dimensions with dimensions="2"
   2. setting the length of box in Z diretion to zero with lz="0"
   3. specifying the position of particles in Z direction as 0.0
   
An example is given::

	mst_version 1.0
		num_particles
			8
		timestep
			0
		dimension
			2
		box
			200.0    200.0    0.0   
		position
			 28.5678528848   -37.9327360252     0.0000000000
			 28.0019705849   -37.1082499897     0.0000000000
			 29.5648198865   -37.8549105956     0.0000000000
			 28.1367681830   -38.8350474902     0.0000000000
			-37.5589154370   -72.8549398355     0.0000000000
			-38.4958248509   -72.5053675968     0.0000000000
			-36.7877222908   -72.2183386015     0.0000000000
			-37.3931991693   -73.8411133084     0.0000000000 
	mst_end	

3. For molgen script to generate a two-dimensional configuration file, a specification of two dimensions 
and box size in Z direction as 0.0 is necessary. Such as::


    import molgen
    
    mol=molgen.Molecule(4)
    mol.setParticleTypes("A,B,B,B")
    mol.setTopology("0-1,0-2,0-3")
    mol.setBondLength("A","B", 1.0)
    mol.setAngleDegree("B", "A", "B", 120)
    mol.setInit("B", 1)
    mol.setCris("A", 1)
    
    
    gen=molgen.Generators(200, 200, 0.0) # box size in X, Y, and Z directions
    gen.addMolecule(mol, 2000)
    gen.setDimension(2)
    gen.setMinimumDistance(1.0)
    gen.outPutMST("pn2d")

  