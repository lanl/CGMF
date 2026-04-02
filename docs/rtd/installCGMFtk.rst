*********************
How to install CGMFtk
*********************

There are two ways of installing `CGMFtk`.  First, pip can be used to install the python package locally.  From the command line within the folder /cgmf/tools/, it can be installed with:

.. code-block:: python

   pip install .

(The pip version used depends on the python version that will be used with the package.  It is recomended that python 3 or above be used.)  

Pip is the recommended way to install `CGMFtk`, as the package will now be installed more globally on the user's machine and can easily be directly imported into a python script or `jupyter_notebook`.

Alternatively (or if pip is not available), the program can be install directly from the setup file:

.. code-block:: python

  python setup.py install 

Again, python 3 and above is recommended.  

This second command will install `CGMFtk` in the /site-packages folder of the corresponding python version directory.  It is important to note that the user does not always have access to local directories, and in this case, the --user flag should be included:

.. code-block:: python

   python setup.py install --user
