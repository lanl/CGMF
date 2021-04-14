================================
CGMFtk - python analysis package
================================

`python` provides a powerful framework to analyze the output of the :program:`CGMF` history files.  We have created a python package, `CGMFtk` that can be used to read the yield file or history file created by :program:`CGMF` and then easily extract results of interest, including calculating neutron and gamma multiplicities, prompt particle energy spectra, isomeric ratios, and correlations between observables of interest.  Details on how to calculate many of these observables are given in the next section.  Here, we show how to install the python package and provide a list of functions within each class.

.. toctree::
    :maxdepth: 1

    Installing CGMFtk<installCGMFtk>
    CGMFtk Classes<classes>
    The History Class<histories>
    The Yield Class<yields>
    