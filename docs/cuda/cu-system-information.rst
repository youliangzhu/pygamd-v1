Information
===========

Perform configuration
---------------------

.. py:class:: PerformConfig(gpu_list)

   The constructor of perform configuration object.
	  
   :param PyObject* gpu_list: The gpu list.

   Example::
   
      global _options
      parser = OptionParser()
      parser.add_option('--gpu', dest='gpu',help='GPU on which to execute')
      (_options, args) = parser.parse_args()
	  
      perform_config = gala.PerformConfig(_options.gpu)



All information
---------------

.. py:class:: AllInfo(reader, perf_conf)
   
   The constructor of all information object.
   
   :param Reader reader: The file parser   
   :param PerformConfig perf_conf: The perform configuration  

   .. py:function:: setNDimensions(unsigned int nds)
   
      set the number of dimensions
	  
   .. py:function:: addParticleType(string ptype)
   
      add a particle type to system
	  
   .. py:function:: addBondType(string bondtype)
   
      add a bond type to system     
      
   .. py:function:: addAngleType(string angletype)
   
      add a angle type to system     
      
   .. py:function:: addBondTypeByPairs()
   
      add bond types to system according to particle types      
      
   .. py:function:: addAngleTypeByPairs()
   
      add angle types to system according to particle types 
   
   Example::
   
      filename = 'PE-1000.xml'
      build_method = gala.XMLReader(filename)
      perform_config = gala.PerformConfig(_options.gpu)
      all_info = gala.AllInfo(build_method, perform_config)


