CGMFtk, CGMF toolkit 
====================

Updated November 2020

---
Versions
--------

### Current Version 1.0

- includes function to analyze CGMF history files and CGMF yield files

---
Installing Instructions
-----------------------

Within the /CGMF/tools/ directory either:
a) Type: pip3 install .
b) Type: python setup.py install
   + the --user flag might have to be included (python setup.py install --user) depending on the permissions of your sysetem

It is recommended that python3 is used both for installation via pip (pip3) or the setup.py file.

From the pip installation, the two classes, CGMFtk.histories() and CGMFtk.yields() should be available anywhere on your system.  If this is not the case (and for the installation from the setup file), before importing the modules use:
   > import sys
   > import site
   > sys.path.append(site.USER_SITE+'/CGMFtk/')

The plotting features in CGMFtkExample.ipynb require python3

---
Execution Instructions
----------------------

The two classes within CGMFtk can be imported into a python script or notebook as
   > from CGMFtk import histories as fh
   > from CGMFtk import yields as yld

A history object is initialized from either class, for example as 
   > hist = fh.Histories('CGMFHistoryFile.cgmf')

Options:
   + file 	[required]	CGMF history or yield file
   + nevents 	[optional]	number of events to be read from the file
		> hist = fh.Histories('CGMFHistoryFile.cgmf',nevents=100)

