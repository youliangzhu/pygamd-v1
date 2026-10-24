Data format
===========

.. _mst-format:

MST format
----------
We take MST format files as the standard input and output configuration files. 
The MST files can contain coordinates, types, masses, velocities, bond connections, angles, dihedrals and so on.
Here is an example of the MST file of a single molecule system. The molecule consisting of four particles is depicted in following picture. 

.. image:: mst-config.png
    :width: 250 px
    :align: center
    :alt: Principle of dihedral torsion

The data in a line of MST file corresponds to a particle and all particles are given in sequence. 
For example, the coordinate of a particle in x, y, and z directions is written in a line and three columns in MST files. 
However, this rule does not include topological relevant information, including bonds, angles and dihedrals.


Snapshot file
^^^^^^^^^^^^^

   An example of MST snapshot file with particles coordinates, velocities, types, masses ... ::

	mst_version 1.0
		num_particles
			4
		timestep
			0
		dimension
			3
		box
			10.0	10.0	10.0
		position
			-1  2 -1
			-2  3  0
			-1  4  1
			-1  5  2
		velocity
			1  2  3
			1  0  0
			3 -2  1
			0  1  1
		type
			A
			B
			B
			A
		mass
			1.0
			2.1
			1.0
			1.0

   The file could include more information, such as bond, angle, dihedral ... :: 
   
     # bond with 'bond type (string), particle index i (int), j (int)'. 
		bond
			polymer 0 1
			polymer 1 2
			polymer 2 3
      
     # angle with 'angle type (string), particle index i (int), j (int), k (int)'. 	  
		angle
			theta 0 1 2
			theta 1 2 3
      
     # dihedral with 'dihedral type (string), particle index i (int), j (int), k (int), l (int)'.
		dihedral
			phi 0 1 2 3

     # virial site with 'vsite type (string), particle index i (int), j (int), k (int), l (int)'.			
		vsite
			v 3 0 1 2			

     # the diameter of particles with float type.
		diameter
			1.0
			1.0
			1.0
			1.0

     # the charge of particles with float type.
		charge
			 1.333
			 1.333
			-1.333
			-1.333

     # the body index of particles with int type, -1 for non-body particles.
		body
			-1
			-1
			 0
			 0
	  
     # the image in x, y, and z directions of particles with int type.	  
		image
			0 0 0 
			0 0 0
			0 0 0
			0 0 0
	  
     # the velocity in x, y, and z directions of particles with float type. 
		velocity
			 3.768     -2.595    -1.874
			-3.988     -1.148     2.800
			 1.570      1.015    -3.167
			 2.441     -1.859    -1.039


     # the orientation vector (x, y, z) of particles with float type.
		orientation
			-0.922     0.085     0.376
			-0.411    -0.637    -0.651
			 0.293     0.892    -0.342
			-0.223     0.084     0.970  

     # the quaternion vector (x, y, z, w) of particles with float type. 	  
		quaternion
			 0.369   0.817  -0.143   0.418
			-0.516  -0.552   0.653   0.024
			-0.521  -0.002   0.131   0.843
			-0.640   0.159  -0.048  -0.749

     # the angular velocity of rotation in x, y, and z directions of particles with float type.	  
		rotation
			-0.640    0.571   -0.512
			-0.744    0.346    0.569
			 0.620   -0.086    0.779
			-0.542    0.319   -0.776	  

    # the moment of inertia in x, y, and z directions of particles with float type.	  
		inert
			1.0 1.0 3.0
			1.0 1.0 3.0
			1.0 1.0 3.0
			1.0 1.0 3.0
			
    # the rotated angles of in x, y, and z directions of particles with float type.	  
		rotangle
			9.478    -1.677    8.239
			8.908    -1.214    8.086
			9.011    -0.653    7.600
			8.993    -0.488    8.331	

    # the initiator indication of particles with int type, 1 for initiator.	  
		init
			0
			1
			0
			1

    # the crosslinking number of particles with int type, 0 for reactable monomer.	  
		cris
			0
			0
			0
			0

    # the molecule index of particles with int type, -1 for free particles.  
		molecule
			0
			0
			1
			1	 	  

   The attribute of anisotropic particles ... ::

    # the particle patch attribute with 'particle type (string), patch number (int)' 
    # followd by 'patch type(string), patch size (float), 
    # patch position vector in x, y, z directions (float)'.
		patch
			B 2
			p1 60  0    0    1
			p1 60  0    0   -1
	  
    # the patch-patch interaction parameter with 'patch type (string), patch type (string), 
    # gamma_epsilon (float), alpha (float)'.	  
		patch_param
			p1 p1 88.0 0.5
	  
    # the particle shape attribute with 'particle type(string), diameter a, diameter b, diameter c, 
    # epsion a, epsion b, epsion c (float)'. The a, b, c are along x, y, z directions in body frame, 
    # respectively.	  
		asphere
			A 1.0 1.0 1.0 3.0 3.0 3.0      
			B 1.0 1.0 3.0 1.0 1.0 0.2
			
    # the end of file.
	mst_end
	
Trajectory file
^^^^^^^^^^^^^^^

   A MST trajectory file could contain multiple frames. The properties in trajectory file are divied into 
   two classes, i.e. invariant data and variant data. The invarant data is only output once, whereas the variant data is output every frame.

   An example of MST trajectory file::
   
	mst_version 1.0
	invariant_data
		num_particles
			4
		dimension
			3
		box
			10.0	10.00	10.0
		bond
			polymer 0 1
			polymer 1 2
			polymer 2 3	  
		angle
			theta 0 1 2
			theta 1 2 3  
		dihedral
			phi 0 1 2 3
		type
			A
			B
			B
			A	
	variant_data
	frame	0
		timestep
			0
		position
			0	0	0
			1	0	0
			2	0	0
			3	0	0
		image
			0	0	0
			0	0	0
			0	0	0
			0	0	0
	frame_end
	frame	1
		timestep
			10000
		position
			0	1	0
			1	1	0
			2	1	0
			3	1	0
		image
			0	0	0
			0	0	0
			0	0	0
			0	0	0
	frame_end			
	frame	2
		timestep
			20000
		position
			0	2	0
			1	2	0
			2	2	0
			3	2	0
		image
			0	0	0
			0	0	0
			0	0	0
			0	0	0   
	frame_end
	