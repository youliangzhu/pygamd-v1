Import libraries
================

import cuda-lib
---------------

Import the pygamd module from cuda libray for the simulations on NVIDIA GPU.

   Examples::

      import cu_gala as gala

	
.. note::
	
       if the cu_gala could not be found by Python, ``import pygamd`` before ``import cu_gala`` may provide some help to find ``cu_gala`` by searching the built-in paths.
       This method is also adapt to some other libraries, such as ``import molgen`` and ``import hip_gala``.

   Examples::

      import pygamd
      import cu_gala as gala


import hip-lib
--------------


Import the pygamd module from hip libray for the simulations on DCU.

   Examples::

      import hip_gala as gala

