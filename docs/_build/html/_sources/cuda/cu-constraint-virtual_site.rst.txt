Virtual site
============

Description:

	GROMACS method

.. py:class:: Vsite(all_info)

   Constructor of a virtual site object.
 
   :param AllInfo all_info: System information.

   .. py:function:: setParams(string type, float a, float b, float c, VST vst)
   
      specifies the virtual site parameters with a, b, c, and virtual site type. The candidates of VST are 'v2', 'v3', 'v3fd', 'v3fad', 'v3out', and 'v4fdn'.	  

   Example::
   
      vs = gala.Vsite(all_info)#virtual interaction sites
      vs.setParams('v', 0.128012065, 0.128012065, 0.0, gala.VST.v3 )
      app.add(vs)



