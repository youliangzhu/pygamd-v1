Data output
===========

Common functions for Data output
--------------------------------

   .. py:function:: setPrecision(int npre)
   
      Set the number of places after decimal point.
	  
   .. py:function:: setHead(int nhead) 
   
      Set the number of places before decimal point.
	  
	  
Collective information
----------------------

.. py:class:: DumpInfo(all_info, comp_info, filename)

   Constructor of an information dump object for a group of particles.
   
   :param AllInfo all_info: System information.	
   :param ComputeInfo comp_info: Object for calculating collective information.	   
   :param str filename: Output file name.  

   .. py:function:: dumpAnisotropy()
   
      Outputs information related to anisotropic particles.
	  
   .. py:function:: dumpVirial(Force object)
   
      Outputs virials of the Force object.
	  
   .. py:function:: dumpPotential(Force object)
   
      Outputs potentials of the Force object.
	  
   .. py:function:: dumpVirialMatrix(Force object)
   
      Outputs virial matrixes including 'virial_xx', 'virial_xy', 'virial_xz', 'virial_yy', 'virial_yz', and 'virial_zz' of the Force object.
	  
   .. py:function:: dumpPressTensor()
   
      Outputs press tensors including 'press_xx', 'press_xy', 'press_xz', 'press_yy', 'press_yz', and 'press_zz' of the system.   
   
   .. py:function:: dumpTypeTemp(string type)
   
      Outputs temperatures of a type of particles.
	  
   .. py:function:: dumpParticleForce(int tag_i)
   
      Outputs the forces of a particle indicated by tag_i in order to trace it's force condition.

   .. py:function:: dumpParticlePosition(int tag_i)
   
      Outputs the positions of a particle indicated by tag_i in order to trace it's position.
	  
   .. py:function:: dumpBoxSize()  
	  
      Outputs box sizes including the box lengths in 'X', 'Y', and 'Z' directions and volume.
	  
   .. py:function:: setPeriod(int period)  
	  
      Period to output data.	  
	  
   Example::
   
      dinfo = gala.DumpInfo(all_info, comp_info, 'data.log')
      dinfo.setPeriod(200)
      app.add(dinfo)

MOL2 dump
---------

.. py:class:: MOL2Dump(all_info, filename)

   Constructor of an object to dump mol2 files.
	  
   :param AllInfo all_info: System information.   
   :param str filename: Output file base name. 	  

   .. py:function:: setChangeFreeType(string type)
   
      specifies the type of free particles which will be changed to be 'F' in output file.
	  
   .. py:function:: deleteBoundaryBond(bool switch)
   
      switches on the function of screening the bonds across the box with 'True'.
	  
   Example::
   
      mol2 = gala.MOL2Dump(all_info, 'particles')
      mol2.setPeriod(100000)
      mol2.deleteBoundaryBond(True)
      app.add(mol2)

XML dump
--------

.. py:class:: XMLDump(all_info, filename)

   Constructor of an object to dump XML files.

   :param AllInfo all_info: System information.   
   :param str filename: Output file base name.

	  
.. py:class:: XMLDump(all_info, group, filename)

   Constructor of an object to dump XML files for a group of particles.

   :param AllInfo all_info: System information.
   :param ParticleSet group: A group of particles.	
   :param str filename: Output file base name.

   .. py:function:: setOutput(PyObject* out_put_list)
   
      indicates the output data type with the candidates::
	  
		['position', 'type', 'velocity', 'mass', 'image', 'force',
		'potential', 'virial', 'virial_matrix', 'charge', 'diameter',
		'body', 'orientation', 'quaternion', 'rotation', 'rotangle',
		'torque', 'inert', 'init', 'cris', 'molecule', 'bond', 'angle',
		'dihedral', 'constraint', 'vsite']
	
	Each data type also could be outputed by a single function as following.
	  
	  
   .. py:function:: setOutputPosition(bool switch)
   
      Outputs particle position (default value is true).
	  
   .. py:function:: setOutputType (bool switch)
   
      Outputs particle type (default value is true).
	  
   .. py:function:: setOutputImage(bool switch)
   
      Outputs particle image.
	  
   .. py:function:: setOutputVelocity(bool switch)
   
      Outputs particle velocity.
	  
   .. py:function:: setOutputMass(bool switch)
   
      Outputs particle mass.
	  
   .. py:function:: setOutputCharge(bool switch)
   
      Outputs particle charge.
	  
   .. py:function:: setOutputDiameter(bool switch)
   
      Outputs particle diameter.
	  
   .. py:function:: setOutputBody(bool switch)
   
      Outputs particle body.
	  
   .. py:function:: setOutputVirial(bool switch)
   
      Outputs particle virial.
	  
   .. py:function:: setOutputForce(bool switch)
   
      Outputs particle force (x, y, z).
	  
   .. py:function:: setOutputPotential(bool switch)
   
      Outputs particle potential.
	  
   .. py:function:: setOutputOrientation(bool switch)
   
      Outputs particle orientation.
	  
   .. py:function:: setOutputQuaternion(bool switch)
   
      Outputs particle quaternion.
	  
   .. py:function:: setOutputRotation(bool switch)
   
      Outputs particle rotation velocity.
	  
   .. py:function:: setOutputRotangle(bool switch)
   
      Outputs particle accumulated rotated angle at the direction of (0, 0, 1) in 3D or (0, 1, 0) in 2D.
	  
   .. py:function:: setOutputTorque(bool switch)
   
      Outputs particle torque.
	  
   .. py:function:: setOutputInert(bool switch)
   
      Outputs particle inert tensor.
	  
   .. py:function:: setOutputInit(bool switch)
   
      Outputs particle initiator indicator.
	  
   .. py:function:: setOutputCris(bool switch)
   
      Outputs particle cross-linking indicator.
	  
   .. py:function:: setOutputBond(bool switch)
   
      Outputs bonds.
	  
   .. py:function:: setOutputAngle(bool switch)
   
      Outputs angles.
	  
   .. py:function:: setOutputDihedral(bool switch)
   
      Outputs dihedrals.

   .. py:function:: setOutputConstraint(bool switch)
 
      Outputs bond constraints. 
		
   .. py:function:: setOutputVsite(bool switch)
   
      Outputs virual sites.
		
   .. py:function:: setOutputLocalForce(Force object)
   
      Outputs particle force(x, y, z and w where w is potential) of a Force object.

   .. py:function:: setOutputLocalVirial(Force object)
   
      Outputs particle virial of a Force object.   
		
   .. py:function:: setOutputLocalVirialMatrix(Force object)
   
       Outputs particle virial matrix of a Force object.  
		
   .. py:function:: setOutputPatch(AniForce object)	
   
      outputs patch information for display in OVITO.   

   .. py:function:: setOutputEllipsoid(BondForceHarmonicEllipsoid object)  
	  
      outputs ellipsoid bond information for display in OVITO.	  
	  
   .. py:function:: setOutputEllipsoid(PBGBForce object)
   
      outputs ellipsoid information for display in OVITO.

   .. py:function:: setOutputEllipsoid(GBForce object)
   
      outputs ellipsoid information for display in OVITO.	  
	  
   
   Example::
   
      xml = gala.XMLDump(all_info, 'particles')
      xml.setOutput(['image', 'bond'])
      xml.setPeriod(100000)
      app.add(xml)

DCD trajectory dump
-------------------

.. py:class:: DCDDump(all_info, filename, overwrite)

   The constructor of a dump object of DCD file.
	  
   :param AllInfo all_info: The system information.	
   :param str filename: The output file name.
   :param bool overwrite: If overwrite the existed DCD file.   
	  
.. py:class:: DCDDump(all_info, group, filename, overwrite)

   The constructor of a dump object of DCD file for a group of particles.
	  
   :param AllInfo all_info: The system information.
   :param ParticleSet group: The group of particles.	
   :param str filename: The output file name. 
   :param bool overwrite: If overwrite the existed DCD file.    

   .. py:function:: unpbc(bool switch)
   
      Outputs particle positions without the application of periodic boundary condition (PBC). Default value is False, i.e. with PBC condition.
	  
   .. py:function:: unwrap(bool switch)
   
      Unwraps the molecules across box boundary due to PBC condition. Default value is False, i.e. wrapping molecules.
	  
   Example::
   
      dcd = gala.DCDDump(all_info, 'particles', True)
      dcd.unwrap(True)
      dcd.setPeriod(100000)
      app.add(dcd)

Binary dump
--------------------

.. py:class:: BinaryDump(all_info, filename)

   Constructor of an object to dump binary files.
	  
   :param AllInfo all_info: System information.
   :param str filename: Output file base name. 	  

   .. py:function:: setOutput(PyObject* out_put_list)
   
      indicates the output data type with the candidates::
	  
		['position', 'type', 'velocity', 'mass', 'image', 'force',
		'potential', 'virial', 'virial_matrix', 'charge', 'diameter',
		'body', 'orientation', 'quaternion', 'rotation', 'rotangle',
		'torque', 'inert', 'init', 'cris', 'molecule', 'bond', 'angle',
		'dihedral', 'constraint', 'vsite']

   .. py:function:: setOutputAll()
   
      Outputs all data.
	  
   .. py:function:: setOutputForRestart()
   
      Outputs data needed for restarting.
	  
   .. py:function:: enableCompression(bool switch)
   
      Compresses output file.
	  
   Example::
   
      binary = gala.BinaryDump(all_info, 'particle')
      binary.setOutput(['image', 'bond'])
      binary.setPeriod(10000)
      app.add(binary)


