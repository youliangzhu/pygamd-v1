Data format
===========

.. _xml-format:

XML format
----------
We take XML format files as the standard input and output configuration files. 
The XML files can contain coordinates, types, masses, velocities, bond connections, angles, dihedrals and so on.
Here is an example of the XML file of a single molecule system. The molecule consisting of four particles is depicted in following picture. 

.. image:: xml-config.png
    :width: 250 px
    :align: center
    :alt: Principle of dihedral torsion

The data in a line of XML file corresponds to a particle and all particles are given in sequence. 
For example, the coordinate of a particle in x, y, and z directions is written in a line and three columns in XML files. 
However, this rule does not include topological relevant information, including bonds, angles and dihedrals.

   An example XML file with particles coordinates, velocities, types, masses ... ::

      <?xml version="1.0" encoding="UTF-8"?>
      <galamost_xml version="1.3">
      <configuration time_step="0" dimensions="3" natoms="4" >
      <box lx="10" ly="10" lz="10"/>
      <position num="4">
      -1 2 -1
      -2 3 0
      -1 4 1
      -1 5 2
      </position>
      <velocity num="4">
      1 2 3
      1 0 0
      3 -2 1
      0 1 1
      </velocity>
      <type num="4">
      A
      B
      B
      A
      </type>
      <mass num="4">
      1.0
      2.1
      1.0
      1.0
      </mass>
      </configuration>
      </galamost_xml>

   The file could include the nodes of bond, angle, dihedral ... :: 
   
     # bond with 'bond type, the index of particle i, j'. 
      <bond num="3">
      polymer 0 1
      polymer 1 2
      polymer 2 3
      </bond>
      
     # angle with 'angle type, the index of particle i, j, k'. 	  
      <angle num="2">
      theta 0 1 2
      theta 1 2 3
      </angle>
      
     # dihedral with 'dihedral type, the index of particle i, j, k, l'. 	  
      <dihedral num="1">
      phi 0 1 2 3
      </dihedral>
         
   The other nodes of XML ... :: 
      
     # the diameter of particles in float type.
      <diameter num="4">
      1.0
      1.0
      1.0
      1.0
      </diameter>

     # the charge of particles in float type.
      <charge num="4">
       1.333
       1.333
      -1.333
      -1.333
      </charge>

     # the body index of particles in int type, -1 for non-body particles.
      <body num="4">
      -1
      -1
      0
      0
      </body>
	  
     # the image in x, y, and z directions of particles in int3 type.	  
      <image num="4">
      0 0 0 
      0 0 0
      0 0 0
      0 0 0
      </image>
	  
     # the velocity in x, y, and z directions of particles in float3 type. 
      <velocity num="4">
       3.768     -2.595    -1.874
      -3.988     -1.148     2.800
       1.570      1.015    -3.167
       2.441     -1.859    -1.039
      </velocity>


     # the orientation vector (x, y, z) of particles in float3 type.
      <orientation num="4">
       -0.922     0.085     0.376
       -0.411    -0.637    -0.651
        0.293     0.892    -0.342
       -0.223     0.084     0.970  
      </orientation>

     # the quaternion vector (x, y, z, w) of particles in float4 type. 	  
      <quaternion num="4">
       0.369   0.817  -0.143   0.418
      -0.516  -0.552   0.653   0.024
      -0.521  -0.002   0.131   0.843
      -0.640   0.159  -0.048  -0.749  
      </quaternion>

     # the angular velocity of rotation in x, y, and z directions of particles in float3 type.	  
      <rotation num="4">
       -0.640    0.571   -0.512
       -0.744    0.346    0.569
        0.620   -0.086    0.779
       -0.542    0.319   -0.776	  
      </rotation>	  

    # the moment of inertia in x, y, and z directions of particles in float3 type.	  
      <inert num="4">
      1.0 1.0 3.0
      1.0 1.0 3.0
      1.0 1.0 3.0
      1.0 1.0 3.0	  
      </inert>	  

    # the initiator indication of particles in int type, 1 for initiator.	  
      <h_init num="4">
      0
      1
      0
      1
      </h_init>	 

    # the crosslinking number of particles in int type, 0 for reactable monomer.	  
      <h_cris num="4">
      0
      0
      0
      0
      </h_cris>	 

    # the molecule index of particles in int type.	  
      <molecule num="4">
      0
      0
      1
      1
      </molecule>	 	  

   The nodes of anisotropic particle attribute ... ::

    # the particle patch attribute with 'particle type, patch number' 
    # followd by 'patch type, patch size, patch position vector in x, y, z directions'.
      <Patches>
      B 2
      p1 60  0    0    1
      p1 60  0    0   -1
      </Patches>
	  
    # the patch-patch interaction parameter with 'patch type, patch type, gamma_epsilon, alpha'.	  
      <PatchParams>
      p1 p1 88.0 0.5
      </PatchParams>
	  
    # the particle shape attribute with 'particle type, diameter a, diameter b, diameter c, 
    # epsion a, epsion b, epsion c'. The a, b, c are along x, y, z directions in body frame, 
    # respectively.	  
      <Aspheres>
      A 1.0 1.0 1.0 3.0 3.0 3.0      
      B 1.0 1.0 3.0 1.0 1.0 0.2 
      </Aspheres>



   
   
   