Particle set
============

.. py:class:: ParticleSet(all_info, keyword)

   The constructor of particle set object with keyword.
  
   :param AllInfo all_info: The system information.
   :param str keyword: Keywords or the python list of keywords. The candidates of keyword are “all”, “body”, “non_body”, “charge”, and particle type(string).

  
.. py:class:: ParticleSet(all_info, min, max)

   The constructor of particle set object with the range of particle index.
  
   :param AllInfo all_info: The system information.
   :param int min: The particle index range from min to max.
   :param int max: The particle index range from min to max.  

.. py:class:: ParticleSet(all_info, members_list)

   The constructor of particle set object with python object.
  
   :param AllInfo all_info: The system information.
   :param PyObject* members_list: The Python object. Note: multiple keywords and particle index can be included in Python object 
      
.. py:class:: ParticleSet(all_info, member_tags)

   The constructor of particle set object with a vector of particle tags.
  
   :param AllInfo all_info: The system information.
   :param vector<unsigned_int> member_tags: A vector of particle tags	  
  
   .. py:function:: ParticleSet combine(ParticleSet group1, ParticleSet group2)
   
      combines two particle groups into one.

   .. py:function:: getNumMembers()
   
      return the number of group members .
   
   Example::
   
      groupC = gala.ParticleSet(all_info,'C')
      # initializes a particle set object by a particle type.
	  
      group = gala.ParticleSet(all_info,['A', 'B', 'C'])
      # initializes a particle set object by particle types	  
	  
      groupB = gala.ParticleSet(all_info,'body')
      # initializes a particle set object of body particles by 'body'.
	  
      group_e = gala.ParticleSet(all_info,'charge')
      # initializes a particle set object of charged particles by 'charge'.
	  
      groupC = gala.ParticleSet(all_info,'all')
      # initializes a particle set object of all particles by 'all'.
	  
      group= gala.ParticleSet(all_info, ['A', 12, 'body'])
      # initializes a particle set object of particles by mixed keywords.
	  
      groupAB= gala.ParticleSet.combine(groupA, groupB)
      # combines two particle sets into one set.
   
   
Dynamic Particle Set
====================

.. py:class:: DynamicParticleSet(all_info, keyword)

   The constructor of particle set object with keyword. The member of the particle set is automatically updated every time step.
  
   :param AllInfo all_info: The system information.
   :param str keyword: Keywords or the python list of keywords. The candidates of keyword are particle type(string).

  
.. py:class:: DynamicParticleSet(all_info, lx_min, lx_max, ly_min, ly_max, lz_min, lz_max)

   The constructor of particle set object with space range. The member of the particle set is automatically updated every time step.
  
   :param AllInfo all_info: The system information.
   :param int lx_min: The particles in the box with lx >= lx_min.
   :param int lx_max: The particles in the box with lx <  lx_max.
   :param int ly_min: The particles in the box with ly >= ly_min.
   :param int ly_max: The particles in the box with ly <  ly_max. 
   :param int lz_min: The particles in the box with lz >= lz_min.
   :param int lz_max: The particles in the box with lz <  lz_max. 
   
   Example::
   
      groupC = gala.DynamicParticleSet(all_info, -10.0, 10.0, -10.0, 10.0, -2.0, 2.0)
      # initializes a particle set object by spatial range.
   