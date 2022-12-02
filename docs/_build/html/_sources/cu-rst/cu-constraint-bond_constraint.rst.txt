Bond constraint
===============

Description:

	LINCS algorithm


.. py:class:: BondConstraint(all_info)

   Constructor of a bond bonstraint object.
 
   :param AllInfo all_info: System information.

   .. py:function:: setParams(string type, float k, float r0)
   
      specifies the bond constraint parameters with bond type and equilibrium length.

   .. py:function:: setNumIters(int ncycles)
   
      specifies the number of iterations of calcuation.
		
   .. py:function:: setExpansionOrder(int order)

      specifies the spread order.   

   Example::
   
      bondconstraint = gala.BondConstraint(all_info)# bond constraints using LINCS algorithm
      bondconstraint.setParams('oh', 0.09572)#(type, r0)
      bondconstraint.setParams('hh', 0.15139)#(type, r0)
      bondconstraint.setExpansionOrder(4)
      bondconstraint.setNumIters(1)
      app.add(bondconstraint)

	  
	  




