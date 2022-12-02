MD-SCF
======

MDSCF force
-----------

Description:

   Using hybrid particle-field technique to accelerate CGMD simulations (G. Milano and T. Kawakatsu, J. Chem. Phys. 130, 214106, 2009). 
   This method could largely speed up some slowly evolving processes in CGMD simulations, such as microphase separation and self-assembly 
   of polymeric systems.

.. py:class:: MDSCFForce(AllInfo all_info, int nx, int ny, int nz, float comp)

   Constructor of an object of MD-SCF force.
   
   :param AllInfo all_info: System information.
   :param int nx: Number of grid in x direction.  
   :param int ny: Number of grid in y direction.  
   :param int nz: Number of grid in z direction.     
   :param float comp: Compressibility.
   
   .. py:function:: setPeriodScf(int idl2_step, int idl_step)
   
      sets the periods of computing density field and updating density field.

   .. py:function:: setParams(string type1, string type2, float chi)
   
      specifies the MD-SCF interaction parameters by pair types with type1, type2, chi parameter.
	  
   .. py:function:: setNewVersion(bool switch)
   
      switches the function of newly developed method of implementation.
	  
   Example::
   
      scf = gala.MDSCFForce(all_info, 22, 22, 22, 0.100)
      scf.setParams('N', 'N',0.000)
      scf.setParams('N', 'P',-1.500)
      scf.setParams('P', 'P',0.000)
      scf.setPeriodScf(1,300)
      scf.setNewVersion(True)
      # switches to newly developed version.
      app.add(scf)	
      # adds this object to the application.



MDSCF electrostatic force
-------------------------

.. py:class:: PFMEForce(AllInfo all_info, int nx, int ny, int nz, float kappa, float epsilonr)

   Constructor of an object to calculate electrostatic forces in MD-SCF framework.
   
   :param AllInfo all_info: System information.
   :param int nx: Number of grid in x direction.  
   :param int ny: Number of grid in y direction.  
   :param int nz: Number of grid in z direction.     
   :param float kappa: kappa = 1/(sqrt(2)sigma) where sigma is the standard deviation of the Gaussian charge distribution.
   :param float epsilonr: Relative dielectric constant.
   
   .. py:function:: setPeriodPFME(int idl2_step, int idl_step)
   
      sets the periods of computing density field and updating density field.
	  
   .. py:function:: setNewVersion(bool switch)
   
      switches the function of newly developed method of implementation.
	  
   Example::
   
      pfme = gala.PFMEForce(all_info,    32,   32,   36,  2.45,  epsilonr)#(mx,my,mz, kappa, epsilonr)
      pfme.setPeriodPFME(1,  100)#(idl2_step, idl_step)
      app.add(pfme)