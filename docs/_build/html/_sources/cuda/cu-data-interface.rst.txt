Data interface
==============
 
System information
------------------

.. py:class:: AllInfo(reader, perform_config)

   Constructor of an all information object for system.
   
   :param Reader reader: reading system information from a file.
   :param PerformConfig perform_config: setting configuration for the simulation.

   .. py:function:: getBasicInfo()
   
      Return an object of basic information.
	  
   .. py:function:: getBondInfo()
   
      Return an object of bond information.
	  
   .. py:function:: getAngleInfo()
   
      Return an object of angle information.
	  
   .. py:function:: getDihedralInfo()
   
      Return an object of dihedral information.
	  
   .. py:function:: getConstraintInfo()
   
      Return an object of constraint information.
	  
   .. py:function:: getVsiteInfo()
   
      Return an object of virtual site information.
	  
   .. py:function:: getPairInfo()
   
      Return an object of pair information.
	  
   .. py:function:: getIntegInfo()
   
      Return an object of integration information.
	  
   .. py:function:: getMoleculeInfo()
   
      Return an object of molecule information.
	  
   .. py:function:: getRigidInfo()
   
      Return an object of rigid body information.
  
   Example::
   
      filename = 'rsCG-modified.xml'
      build_method = gala.XMLReader(filename)
      perform_config = gala.PerformConfig(_options.gpu)
      all_info = gala.AllInfo(build_method, perform_config)

Python interface of BasicInfo
-----------------------------

.. py:class:: BasicInfo(all_info)

   Constructor of a basic information object of system.
	  
   :param AllInfo all_info: System information.
   
   .. py:function:: getParticleTypeMap();

      return a python list of type map
	  
   .. py:function:: getParticlePos(); 
   
      return a python numpy array of position with (N, 4)

   .. py:function:: getParticleForce()
   
      return a python numpy array of Force with (N, 4)
	  
   .. py:function:: getParticleCharge()
   
      return a python numpy array of Charge with (N, 4)
	  
   .. py:function:: getParticlePos_list()
   
      return a python list of position with (N, 4)
	  
   .. py:function:: getParticleForce_list()
   
      return a python list of Force with (N, 4)
	  
   .. py:function:: getParticleCharge_list()
   
      return a python list of Charge with (N, 4)
	  
   .. py:function:: getParticleInit()
   
      return an numpy array of Init with (N, 1)
	  
   .. py:function:: getParticleCris()
   
      return an numpy array of Cris with (N, 1)
	  
   .. py:function:: addParticleType(string name)
   
      add a new particle type to the system

   .. py:function:: setParticleType(unsigned int tag, string name)
   
      set the type of a particle with tag index

   .. py:function:: setParticlePos(unsigned int tag, float px, float py, float pz)
   
       set the position (px, py, pz) of a particle with tag index

   .. py:function:: setParticleForce(unsigned int tag, float fx, float fy, float fz)
   
       set the force (fx, fy, fz) of a particle with tag index

   .. py:function:: setParticleCharge(unsigned int tag, float ci)

       set the charge (ci) of a particle with tag index

   .. py:function:: setParticleInit(unsigned int tag, unsigned int ii)
   
        set the init (ii) of a particle with tag index

   .. py:function:: setParticleCris(unsigned int tag, unsigned int ci)
   
        set the cris (ci) of a particle with tag index

   Example::
   
      all_info.initBasicInfo()
      bi = all_info.getBasicInfo()
      pos = bi.getParticlePos()
      print(pos[0])
      bi.setParticlePos(1, 0.1, 0.2, 0.4)
      pos = bi.getParticlePos()
      print(pos[0])