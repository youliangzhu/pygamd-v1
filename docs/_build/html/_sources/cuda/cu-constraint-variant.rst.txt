Variant
=======

Variant Const
-------------
Description:

    This variant returns a constant value.

.. py:class:: VariantConst(value)

   The constructor of a constant value method.
   
   :param float value: The constant value.
	  
   Example::
   
      v = gala.VariantConst(1.0)
      # set the constant value.

Variant Linear
--------------

    This variant returns a linear interpolation between v0 and v1, or v1 and v2, or until vn-1 and vn

.. py:class:: VariantLinear()

   The constructor of a linearly varying value method.

  .. py:function:: setPoint(unsigned int timestep, double value)
  
      specifies the value at the time step.
	  
   Example::
   
      v = gala.VariantLinear()
      v.setPoint(0, 1.0)         # time step, v0
      v.setPoint(100000, 2.0)    # time step, v1
      v.setPoint(200000, 1.0)    # time step, v2
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

Variant Rsqrt
-------------
Description:

    This variant is usually used for the tensile test at constant volume. When the one side of the length of box is linearly changed from v0 to v1, the changes of the other two sides of box can be calculated by this variant.
    This method return factor*sqrt(v0/v) , where v is real time value of linear interpolation between v0 and v1.

.. py:class:: VariantRsqrt()

   The constructor of a reversed sqrt() varying object.

  .. py:function:: setPoint(unsigned int timestep, double value)
  
   specifies the value at the time step.
   
  .. py:function:: setFactor(double factor)
  
   specifies the factor of sqrt.
   
   Example::
   
      # The box size is (lx = 30, ly = 60, lz = 120) and the lz is stretched from 120 to 240. 
      # To keep the constant volume of box, the changes of lx and ly can be calculated as following.
      
      vz = gala.VariantLinear()  # the change of lz
      vz.setPoint(0, 120)        #time step, box length.
      vz.setPoint(80000000, 240)
      
      vx = gala.VariantRsqrt()   # the change of lx
      vx.setPoint(0, 120)        # time step, v0
      vx.setPoint(80000000, 240) # time step, v1
      vx.setFactor(30)
      
      v2 = gala.VariantRsqrt()   # the change of ly
      v2.setPoint(0, 120)        # time step, v0
      v2.setPoint(80000000, 240) # time step, v1
      v2.setFactor(60)

      # The stretching method to control the changes of box at three directions
      axs = gala.AxialStretching(all_info, group)
      axs.setBoxLength(vz, 'Z')
      axs.setBoxLength(vx, 'X')
      axs.setBoxLength(vy, 'Y')
      axs.setPeriod(1000)
      app.add(axs)

