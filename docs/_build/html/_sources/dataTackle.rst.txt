dataTackle
==========

Usage
-----

The **dataTackle** is a plugin of PYGAMD to analyze some important properties by reading the generated configuration files. 
The **dataTackle** can be compiled and installed by "sh compile.sh". You can use this tool to analyze one or more files at a time.

   Examples for running::
   
      dataTackle particle.mst
      dataTackle particle0.mst particle1.mst
      dataTackle *.mst
	  
After pressing enter, a menu of function options will be listed. You can choose one or more functions by 
the number indexes separated by blank. The parameters for a function can be input after ":" and seperated by "|". 

   Such as::
   
      4:npot=1000 3:gpu=0|rmax=2.0

The **dataTackle** also can be called and passed parameters in a Shell script.

   Such as::
   
      echo -e "18:gpu=1|rmax=3.0" | dataTackle particles.*.mst

The **dataTackle** plugin supports the configuration files with MST and DCD formats. 
The MST files can be tackled independently. However, a DCD trajectory file has to be tackled along with 
a MST file for particles attributes and topological information.

   Examples::

      dataTackle particle.mst trajectory.dcd

For the help about the parameters, you could input "function number:h" after the function list shown.
 	  
   Examples::

      14:h	  
	  
Functions
---------

   Function list::
   
      -------------------------------------------------------------
      1  Rg^2              2  Ed^2               3  RDF              
      4  bond_distri       5  angle_distri       6  dihedral_distri  
      7  stress tensor     8  density            9  unwrapping       
      10 MSD               11 RDF-CM             12 MSD-CM           
      13 ents              14 strfac             15 domain size      
      16 dynamic strfac    17 config check       18 RDF between types
      19 MST conversion    20 patch/spot display 21 SSF               
      22 ADF               23 CND                24 MSAD              
      25 RMSAD             26 ISF                27 OACF              
      28 Q4Q6              29 VORONOI            30 NGP               
      31 RNGP              32 VHF                33 RVHF              
      34 fpSus             35 RfpSus             36 OvlaF             
      37 CISF              38 CAGEISF            39 CAGEMSD           
      40 RMSD              41 P2P4               42 CRYSTALLINITY     
      43 G6_3D             44 W4W6
      -------------------------------------------------------------

1  Rg^2:
^^^^^^^^

   Description:
      The mean square of gyration radius is calculated and output to ``rg2.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		R{_{g}^{2}}&=&\frac{1}{N}\sum\limits_{i=1}^{N}{{{({{{\vec{R}}}_{i}}-{{{\vec{R}}}_{cm}})}^{2}}} \\
		{{\vec{R}}_{cm}}&=&\frac{1}{N}\sum\limits_{i=1}^{N}{{{{\vec{R}}}_{i}}}
        \end{eqnarray*}

    Coefficients:

    - :math:`{\vec{R}}_{i}` - monomer position vector 

2  Ed^2:	  
^^^^^^^^
   
   Description:
      The mean square of end-to-end distance is calculated and output to ``ed2.log``.
	  
    .. math::
        :nowrap:

        \begin{eqnarray*}
		E{_{d}^{2}}={( {{{\vec{R}}}_{0}}-{{{\vec{R}}}_{N-1}} )}^{2}
        \end{eqnarray*}

    Coefficients:

    - :math:`{\vec{R}}_{i}` - monomer position vector 	  
	  
3  RDF:	  
^^^^^^^
   
   Description:
      The radial distribution function of all particles is calculated and output to ``filename.rdf``.
      Averaged value among files will be output to ``rdf.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false
	  
4  bond_distri:	  
^^^^^^^^^^^^^^^

   Description:
      The distribution of bond lengths is calculated and output to ``bond_distr.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		bond\_distri(i \cdot dr)=N(i)/(N \cdot dr)
        \end{eqnarray*}

    Coefficients:

    - :math:`dr` - the space of bond length `L/(2npot)`, where `L` is the box size
    - :math:`N(i)` - the number of bonds in the range of `idr < r < (i+1)dr`, where `i` is an integer
    - :math:`N` - the total number of bonds		
	  
   Parameters:
      :npot=2001

5  angle_distri:	  
^^^^^^^^^^^^^^^^
   
   Description:
      The distribution of angle degrees is calculated and output to ``angle_distr.log``.
	  
    .. math::
        :nowrap:

        \begin{eqnarray*}
		angle\_distri(i \cdot da)=N(i)/(N \cdot da)
        \end{eqnarray*}

    Coefficients:

    - :math:`da` - the space of angle radian `pi/npot`
    - :math:`N(i)` - the number of angles in the range of `ida < angle < (i+1)da`, where `i` is an integer	  
    - :math:`N` - the total number of angles	
	
   Parameters:
      :npot=2001
	  
6  dihedral_distri:	  
^^^^^^^^^^^^^^^^^^^
   
   Description:
      The distribution of dihedral degrees is calculated and output to ``dihedral_distr.log``.

    .. math::
        :nowrap:

        \begin{eqnarray*}
		dihedral\_distri(i \cdot da)=N(i)/(N \cdot da)
        \end{eqnarray*}

    Coefficients:

    - :math:`da` - the space of dihedral angle radian `2pi/npot`
    - :math:`N(i)` - the number of dihedrals in the range of `ida < dihedral angle < (i+1)da`, where `i` is an integer
    - :math:`N` - the total number of dihedrals		
	  
   Parameters:
      :npot=2001
	  
7  stress tensor:	  
^^^^^^^^^^^^^^^^^
   
   Description:
      Stress tensor is calculated by inputing the parameters for force calculation. Result will be output to ``stress_tensor.log``.
	  
   Parameters:
      :bondex=true|bodyex=true|diameter=true 

8  density:	  
^^^^^^^^^^^
   
   Description:
      Real density (g/cm^3) with basic units [amu] and [nm] is calculated and output to ``density.log``.
	  
9  unwrapping:	  
^^^^^^^^^^^^^^
   
   Description:
      This function would unwrap or shift molecules by changing image or position of particles. New configuration 
      will be output to ``filename.reimage.mst``
	  
   Parameters:
      :unwrap_molecule=true|label_free_particle=particle type|molecule_center_in_box=false|
      shiftx=0.0|shifty=0.0|shiftz=0.0|remove_image=false|add_image_to_pos=true|
      remove_bond_cross_box=false|body_keep=false

	  
10 MSD:	  
^^^^^^^
   
   Description:
      The mean square displacement of all particles is calculated and output to ``msd.log``.

   Parameters:
      :direction=XYZ (candidates are X,Y,Z,XY,YZ,XZ,XYZ) 
	  
11 RDF-CM:	  
^^^^^^^^^^
   
   Description:
      The radial distribution function of the mass center of molecules is calculated and output to ``rdf_cm.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2	  
	  
12 MSD-CM:	  
^^^^^^^^^^
   
   Description:
      The mean square displacement of the mass center of molecules is calculated and output to ``msd_cm.log``.
	  
   Parameters:
      :direction=XYZ (candidates are X,Y,Z,XY,YZ,XZ,XYZ)  
	  
13 ents:	  
^^^^^^^^
   
   Description:
      This function would analyze the entanglements of polymers. Result will be output to ``ents.log``.
	  
14 strfac:	  
^^^^^^^^^^
   
   Description:
      The structure factor of particles is calculated and output to ``filename.strf``.
      The averaged value among files is output to ``strf.log``.
	  
   Parameters:
      :qmax=160pi/Lmin|gpu=0|deltaq=2pi/Lmin|direction=XYZ|2D=false

15 domain size:	  
^^^^^^^^^^^^^^^
   
   Description:
      The domain size of components in mixtures is calculated and output to ``domsize.log``.
	  
   Parameters:
      :kmax=20|qc=0.4690|gpu=0

16 dynamic strfac:	  
^^^^^^^^^^^^^^^^^^
   
   Description:
      Dynamic structure factor (incoherent intermediate) measures the decorrelation of the positions 
      of individual monomers with the time on length scale :math:`1/q`, where :math:`q=2\pi\sqrt{x^2+y^2+z^2}/L`, and :math:`L` is cubic box length. 
      :math:`\mbox{kmax}` limits the space in which the :math:`q` with possible combinations of x, y, z will be generated.

      Results will be output to ``dstrf.log``.

   Parameters:
      :kmax=int(L)|q=7.0
	  
   `Maintainer`: Shu-Jia Li
	  
17 config check:	  
^^^^^^^^^^^^^^^^
   
   Description:
      This function would check the configuration including the minimum distance of particles, and the maximum and minimum length of bonds, etc.
      Result will be output to ``config_check.log``.

   Parameters:
      :bondex=true|angleex=true|dihedralex=true|bodyex=true|rcut=2.0


18 RDF between types:	  
^^^^^^^^^^^^^^^^^^^^^
   
   Description:
      The radial distribution function between types is calculated and output to ``filename.type.rdf``.
      The averaged value among files will be output to ``rdf_by_type.log``.
	  
   Parameters:
      :maxbin=100|gpu=0|rmax=Lx/2|bondex=false|angleex=false|molex=false


19 MST conversion:
^^^^^^^^^^^^^^^^^^
   
   Description:
      This function would convert MST files into new files with another format.
	  
   Parameters:
      :lammps=false|gromacs=false|xml=false
  
	  
 	  