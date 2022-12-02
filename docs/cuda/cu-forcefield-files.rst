Read force parameters
=====================

Force field format
---------------------------

Force field for non-boned interactions::

      <pair_params>
      particle_type1  particle_type2  epsilon  sigma  alpha
      </pair_params>

For bond, angle, and dihedral interactions::
      
      <bond_params>
      bond_type  spring_constant  equilibrium_length  function_type
      </bond_params>
      
      <angle_params>
      angle_type  spring_constant  equilibrium_length  function_type
      </angle_params>
      
      <dihedral_params>
      dihedral_type  parameter1  parameter2 ...  function_type
      </dihedral_params>
	  
For bond constraint and virtual site::

      <constraint_params>
      constraint_type  equilibrium_length  function_type
      </constraint_params>
      
      <vsite_params>
      virtual_site_type  a  b  c  function_type
      </vsite_params>  

An example of force field file::
   
      <pair_params>
      Qa Qa 2.3 0.6 1.0
      Qa Q0 2.0 0.6 1.0
      Qa Na 0.5 0.47 1.0
      Qa C4 0.5 0.47 1.0
      Qa C1 2.0 0.62 1.0
      Qa SC3 0.5 0.47 1.0
      Qa SC1 2.0 0.62 1.0
      Qa SP1c 2.7 0.47 1.0
      Q0 Q0 2.0 0.6 1.0
      Q0 Na 0.5 0.47 1.0
      Q0 C4 0.5 0.47 1.0
      Q0 C1 2.0 0.62 1.0
      Q0 SC3 0.5 0.47 1.0
      Q0 SC1 2.0 0.62 1.0
      Q0 SP1c 2.7 0.47 1.0
      Na Na 2.3 0.47 1.0
      Na C4 2.7 0.47 1.0
      Na C1 2.7 0.47 1.0
      Na SC3 2.7 0.47 1.0
      Na SC1 2.7 0.47 1.0
      Na SP1c 2.3 0.47 1.0
      C4 C4 4.5 0.47 1.0
      C4 C1 4.0 0.47 1.0
      C4 SC3 4.5 0.47 1.0
      C4 SC1 4.0 0.47 1.0
      C4 SP1c 2.7 0.47 1.0
      C1 C1 4.5 0.47 1.0
      C1 SC3 4.5 0.47 1.0
      C1 SC1 4.5 0.47 1.0
      C1 SP1c 2.3 0.47 1.0
      SC3 SC3 3.4 0.43 1.0
      SC3 SC1 3.4 0.43 1.0
      SC3 SP1c 2.7 0.47 1.0
      SC1 SC1 3.4 0.43 1.0
      SC1 SP1c 2.3 0.47 1.0
      SP1c SP1c 2.3 0.47 1.0
      </pair_params>
      
      <constraint_params>
      SP1c-SC3 0.4904 1
      SP1c-SC1 0.6019 1
      SC3-SC1 0.2719 1
      SC1-SC3 0.7237 1
      SC1-SC1 0.5376 1
      </constraint_params>
      
      <vsite_params>
      SC1-SC1-SC3-SC1 0.9613 0.6320 0.0 1
      SC1-SC3-SP1c-SC1 0.5207 0.2882 -1.03168 4
      SC1-SC1-SC3-SC1_1 0.2287 0.4111 1.41920 4
      </vsite_params>
      
      <bond_params>
      Q0-Qa 1250.0 0.450 1
      Qa-Na 1250.0 0.450 1
      Na-Na 1250.0 0.370 1
      Na-C1 1250.0 0.480 1
      C1-C1 1250.0 0.480 1
      C1-C4 1250.0 0.480 1
      C4-C4 1250.0 0.480 1
      C4-C1 1250.0 0.480 1
      SC1-C1 1250.0 0.425 1
      </bond_params>
      
      <angle_params>
      Qa-Na-Na 25.0 120.000 2
      Qa-Na-C1 25.0 180.000 2
      Na-C1-C1 35.0 180.000 2
      C1-C1-C1 35.0 180.000 2
      Na-C1-C4 35.0 180.000 2
      C1-C4-C4 20.0 95.000 2
      C4-C4-C1 45.0 120.000 2
      SC1-SC1-C1 25.0 180.0 2
      </angle_params>
      
      <dihedral_params>
      SP1c-SC3-SC1-SC1_F2 -179.7 50.0 2 
      </dihedral_params>
  
Use force fields
-------------------------

Description:

   Force fields in the format could be read by ``force_field_gala`` module. The classes of ``force_field_gala`` module are listed as following.

.. py:class:: LJCoulombShiftForce(all_info, nlist, rcut, rshift, epsilon_r, file)

   Constructor of an object to simutaneously calculate modified Lennard-Jones and Coulomb interactions which are smoothed by a shift function same to GROMACS.
   
   :param AllInfo all_info: System information.
   :param NeighborList nlist: Neighbor list.  
   :param float rcut: Cut-off radius.
   :param float rshift: Shift radius.   
   :param float epsilon_r: Relative dielectric constant.
   :param string file: Force field file.   

   Example::
   
     import force_field_gala
     e_r = 15.0
     lj = force_field_gala.LJCoulombShiftForce(all_info, nlist, 1.2, 0.9, e_r, "Equ.force_field")
     app.add(lj)
	  
.. py:class:: LJEwaldForce(all_info, nlist, rcut, file)

   Constructor of an object to simutaneously calculate Lennard-Jones and the short-part Coulomb interactions.
   
   :param AllInfo all_info: System information.
   :param NeighborList nlist: Neighbor list.  
   :param float rcut: Cut-off radius.
   :param string file: Force field file.    

   .. py:function:: setEnergy_shift()
   
      calls the function to shift LJ potential to be zero at cut-off point.
	  
   .. py:function:: setDispVirialCorr(bool open)
   
      switches the dispersion virial correction.

   Example::

     import force_field_gala 
     lj = force_field_gala.LJEwaldForce(all_info, neighbor_list, 1.0, "ffnonbonded.force_field")
     lj.setEnergy_shift()
     lj.setDispVirialCorr(True)#dispersion virial correction
     app.add(lj)

.. py:class:: BondForceHarmonic(all_info, file)

   Constructor of an object to calculate harmonic bond interactions.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

     bondforce = force_field_gala.BondForceHarmonic(all_info, "ffbonded.force_field")
     app.add(bondforce)

.. py:class:: AngleForceHarmonicCos(all_info, file)

   Constructor of an object to calculate harmonic cosine angle interactions.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

      angleforce = force_field_gala.AngleForceHarmonicCos(all_info, "ffbonded.force_field")
      app.add(angleforce)

.. py:class:: AngleForceHarmonic(all_info, file)

   Constructor of an object to calculate harmonic angle interactions.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

      angleforce = force_field_gala.AngleForceHarmonic(all_info, "ffbonded.force_field")
      app.add(angleforce)

.. py:class:: DihedralForceAmberCosine(all_info, file)

   Constructor of an object to calculate Amber cosine dihedral interactions.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

      dihedralforce = force_field_gala.DihedralForceAmberCosine(all_info, "ffbonded.force_field")
      app.add(dihedralforce)

.. py:class:: DihedralForceHarmonic(all_info, file)

   Constructor of an object to calculate harmonic dihedral interactions.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

      dihedralforce = force_field_gala.DihedralForceHarmonic(all_info, "ffbonded.force_field")
      app.add(dihedralforce)
	  
.. py:class:: BondConstraint(all_info, file)

   Constructor of an object to implement bond constraints.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   .. py:function:: setNumIters(int ncycles)

      specifies the number of iterations of calcuation.

   .. py:function:: setExpansionOrder(int order)

      specifies the spread order.   

   Example::

      bond_constraint = force_field_gala.BondConstraint(all_info, "Equ.force_field")
      bond_constraint.setExpansionOrder(4)
      bond_constraint.setNumIters(1)
      app.add(bond_constraint)
	  
	  
.. py:class:: Vsite(all_info, file)

   Constructor of an object to implement virtual sites using a same method to GROMACS.

   :param AllInfo all_info: System information.
   :param string file: Force field file.

   Example::

      vs = force_field_gala.Vsite(all_info, "Equ.force_field")
      app.add(vs)
  
Use GROMACS force fields
------------------------

Description:

   Force fields in GROMACS format are supported by ``force_field_itp`` module. The usage and methods are same to ``force_field_gala`` module, but for reading the force fields in GROMACS format from itp files.

An example::

      import force_field_itp
      lj = force_field_itp.LJEwaldForce(all_info, neighbor_list, 1.0, "ffnonbonded.itp")
      lj.setEnergy_shift()
      lj.setDispVirialCorr(True)#dispersion virial correction
      app.add(lj)
	  
Convert GROMACS files
---------------------

Description:

   Convert GROMACS files to GALA files including configuration and force fields by ``gro_to_xml`` module.
   Execution command is ``python gro_to_xml.py`` with two necessary parameters ``--gro=`` and ``--top=`` to set
   the GROMACS file names.

An example::

      python gro_to_xml.py --top=Topol.top --gro=Equ.gro  
	  
Then two files 'Equ.xml' of configuration and 'Equ.force_field' of force field will be generated. 