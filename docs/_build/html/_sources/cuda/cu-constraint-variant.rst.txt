Variant
=======

Variant Const
-------------

.. py:class:: VariantConst(value)

   The constructor of a constant value method.
   
   :param float value: The constant value.
	  
   Example::
   
      v = gala.VariantConst(1.0)
      # set the constant value.

Variant Linear
--------------

.. py:class:: VariantLinear()

   The constructor of a linearly varying value method.

  .. py:function:: setPoint(unsigned int timestep, double value)
  
      specifies the value at the time step.
	  
   Example::
   
      v = gala.VariantLinear()
      v.setPoint(0, 1.0)
      v.setPoint(100000, 2.0)
      # set the value at the time step. The value at a time step 
      # varies by linear interpolation.

Variant Sin
-----------

.. py:class:: VariantSin()

   The constructor of a sinusoidal curve varying object.

  .. py:function:: setPoint(unsigned int timestep, double period, double ubd, double lbd)
  
      Function: specifies the period, upper, and lower bounds at the time step.
	  
   Example::
   
      v = gala.VariantSin()
      v.setPoint(0, 1000, 1.0, -1.0)
      v.setPoint(100000, 1000, 2.0, -2.0)
      # set the parameters of sinusoid at the time step and the parameters 
      # at any time step can be gotten by linear interpolation.

Variant Well
------------

.. py:class:: VariantWell()

   The constructor of a well curve varying object.

  .. py:function:: setPoint(unsigned int timestep, double period, double ubd, double lbd)
  
   specifies the period, upper, and lower bounds at the time step.
   
   Example::
   
      v = gala.VariantWell()
      v.setPoint(0, 1000, 1.0, -1.0)
      v.setPoint(100000, 1000, 1.0, -1.0)
      # set the parameters of periodic well at the time step and the parameters 
      # at any time step can be gotten by linear interpolation.


