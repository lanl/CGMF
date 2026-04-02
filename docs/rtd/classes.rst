**********************
Classes within CGMFtk
**********************


The first class is the Histories class, which can be imported with

.. code-block:: python

   from CGMFtk import histories as fh

Note that if `CGMFtk` was installed using the python setup file, you will have to append the location of the package to your python path before importing the `histories` class,

.. code-block:: python
   
   import sys
   sys.path.append('location-of-CGMFtk')

To load a :program:`CGMF` history file into a history object

.. code-block:: python
 
   hist = fh.Histories('history.cgmf')

The second class is the Yields class which can be imported with

.. code-block:: python

   from CGMFtk import yields as yld

Loading a :program:`CGMF` yield file is much the same as loading a history file

.. code-block:: python

   yields = yld.Yields('yields.cgmf')

Note that for both classes, you can provide a relative or an absolute path to the history file, if the file is not in the current directory.
