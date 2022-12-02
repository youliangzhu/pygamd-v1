Application
===========

Modules management
------------------

Usually, we only define an application object in the context of script. The modules can be 
added ``add()`` into or removed ``remove()`` from the application before runing ``run()`` the simulation. 

.. py:class::  Application(all_info, dt)

   The constructor of an application object.
	  
   :param AllInfo all_info: system information
   :param float dt: integration time step	   

   .. py:function:: add(boost::shared_ptr<*> object)
   
      adds an object to the application.
	  
   .. py:function:: remove(boost::shared_ptr<*> object)
   
      removes an added object.
	  
   .. py:function:: clear()
   
      removes all objects from the application.
	  
   .. py:function:: setDt(float dt)
   
      sets integration time step.	  
	  
   .. py:function:: run(unsigned int N)
   
      runs the simulation for N time steps.
	  
   Example::
   
      dt = 0.001
      app = gala.Application(all_info, dt)
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

   dt = 0.001
   app = gala.Application(all_info, dt)
   app.add(lj)
   app.add(nvt)
   app.add(xml)
   app.run(1000)
  
* Second state simulation::

   app.remove(lj)
   app.remove(nvt)
   app.add(harmonic)
   app.add(npt)
   app.run(1000)
   
   
Two-dimensional simulation 
--------------------------

1. Controlling script, i.e. 'file.gala' script is same for two- and three-dimensional simulations.
 
2. However, configuration file i.e. XML file should indicate two-dimensional 
system by:

   1. pointing out dimensions with dimensions="2"
   2. setting the length of box in Z diretion to zero with lz="0"
   3. specifying the position of particles in Z direction as 0.0
   
An example is given::

   <?xml version="1.0" encoding="UTF-8"?>
   <galamost_xml version="1.3">
   <configuration time_step="0" dimensions="2" natoms="8" >
   <box lx="200" ly="200" lz="0"/>
   <position num="8">
       28.5678528848   -37.9327360252     0.0000000000
       28.0019705849   -37.1082499897     0.0000000000
       29.5648198865   -37.8549105956     0.0000000000
       28.1367681830   -38.8350474902     0.0000000000
      -37.5589154370   -72.8549398355     0.0000000000
      -38.4958248509   -72.5053675968     0.0000000000
      -36.7877222908   -72.2183386015     0.0000000000
      -37.3931991693   -73.8411133084     0.0000000000 
   </position>
   </configuration>
   </galamost_xml>	  

3. For molgen script to generate a two-dimensional configuration file, a specification of two dimensions 
and box size in Z direction as 0.0 is necessary. Such as::


    #!/usr/bin/python
    import sys
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
    gen.outPutXml("pn2d")

  