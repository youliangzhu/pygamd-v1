Usage
=====

With a prepared script, you could run pygamd MD engine for obtaining trajectory.

   Examples::
   
      python3 yourscript.py --gpu=0 >a.log&
	  
Where you could specify the GPU id (default value is 0) with the ``--gpu=`` option and output the screen information into ``a.log`` file.

Here is an example of script for DPD simulation. 

Firstly, importing the pygamd module installed as a package of python3 and reading system information by :py:class:`snapshot.read` from a mst file 

   Examples::

      import pygamd
      mst = pygamd.snapshot.read("AB.mst")
	  
After that, we need to build up an application by :py:class:`application.dynamics` which will call defined and added objects.

   Examples::
   
      app = pygamd.application.dynamics(info=mst, dt=0.04)

Further, we should define objects by the classes of pygamd and pass them to the application, such as the following example: DPD force :py:class:`force.dpd`
, NVT thermostat with GWVV algorithm :py:class:`integration.gwvv`, and th dump of system collective information :py:class:`dump.data`.
      
   Examples::
  	  
      fn = pygamd.force.dpd(info=mst, rcut=1.0)
      fn.setParams(type_i="A", type_j="A", alpha=25.0, sigma=3.0)
      fn.setParams(type_i="A", type_j="B", alpha=40.0, sigma=3.0)
      fn.setParams(type_i="B", type_j="B", alpha=25.0, sigma=3.0)
      app.add(fn)
      
      
      gw = pygamd.integration.gwvv(info=mst, group='all')
      app.add(gw)
      
      di = pygamd.dump.data(info=mst, group='all', file='data.log', period=500)
      app.add(di)

	  
Finally, running the simulation with the number of time steps.
      
   Examples::
  	  
      app.run(10000)


