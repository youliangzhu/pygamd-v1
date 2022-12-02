Data input
==========

.. py:class:: Reader

   The basic class of :py:class:`XmlReader` and :py:class:`BinaryReader`.


XML reader
----------

.. py:class:: XMLReader(filename)

   The constructor of XML file parser object.
	  
   :param str filename: The input XML file name.

   Example::
   
      filename = 'dppc.xml'
      # sets the name of input XML file.
      build_method = gala.XMLReader(filename)
      # builds up a parser object for the input XML file.
   
Binary reader
-------------

.. py:class:: BinaryReader(filename)

   The constructor of binary file parser object.
	  
   :param str filename: The input binary file name.
  
   Example::
   
      filename = 'initial.bin'
      # sets the name of binary file.
	  
      build_method = gala.BinaryReader(filename)
      # builds up a reading object for input binary file.
   
   
   