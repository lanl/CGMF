***************
Discrete Levels
***************

The nuclear structure information needed to describe each fission fragment is obtained from the `RIPL-3 <https://www-nds.iaea.org/RIPL-3/>`_ library. Nuclear structure data are always evolving, depending on the availability of new nuclear structure experiments. :program:`CGMF` does not calculate the excited states in a given nucleus, but instead fully depends on the `ENSDF <http://www.nndc.bnl.gov/ensdf/>`_ or somewhat equivalenty, RIPL3 libraries. The study of isomeric ratios or more generally specific :math:`\gamma` decay chains strongly depends on the quality of the underlying nuclear structure information. 

Preparing the nuclear structure file for :program:`CGMF`
========================================================


A special discrete level file (``cgmfDiscreteLevels.dat``) for use with :program:`CGMF` is produced from the RIPL-3 library, and includes adjustments for incomplete information. The Jupyter notebook ``transformLevelDataFile.ipynb`` and the python class ``Nucleus.py`` are used for this purpose.

#. For each nucleus, the level data are first read from the RIPL-3 complete datafile;
#. The levels are then "fixed" for missing or unassigned spin and parity;
#. Finally, the gamma transitions are "fixed".

Fixing the level scheme
-----------------------

If the spin or/and parity of a level is negative, then the level is kept and its (:math:`J,\pi`) values are chosen according to the distribution assumed in the continuum following the Kawano-Chiba-Koura (KCK) level density systematics [`J. Nucl. Sci. and Tech. 43, 1 (2006) <https://www.tandfonline.com/doi/abs/10.1080/18811248.2006.9711062>`_]. Assuming that the parities are evenly distributed, the missing level parity is chosen randomly.

In some cases, RIPL-3 would provide choices of spin and parity for a uncertain level. In that case, the first option will be selected. 

In the case of the ground-state being completely unknown, default values will be selected according to:

* even-even: :math:`0^+`

* odd-odd: :math:`1^+`

* even-odd: :math:`1/2^+`

.. warning:: Such arbitrary decisions can have an non-negligible impact when studying specific :math:`\gamma` decay chains in a particular nucleus. In that case, extra caution should be put in interpreting the results of :program:`CGMF` by studying the origin of the nuclear structure information available (or not) for this nucleus.


Fixing the :math:`\gamma` transitions
-------------------------------------

In addition to incomplete level schemes, the nuclear structure data from RIPL-3/ENSDF can have incomplete decay chains. 

In the case of the first excited state, if no decay data is available, it is assumed that the level decays 100% into the ground-state emitting a :math:`\gamma` ray of energy corresponding to the excitation energy of this first level.

.. todo:: is this the right thing to do? Not always!

In the case of higher excited levels, the spin and parity are chosen according to the (:math:`J,\pi`) values of levels below in energy, and to which the current uncomplete level could be reached through an E1 transition. 

.. note:: More sophisticated level decay schemes should be employed.


References
==========

#. T.Kawano, S.Chiba and H.Koura, *Phenomenological Nuclear Level Densities using the KTUY05 Nuclear Mass Formula for Applications Off-Stability*, `J. Nucl. Sci. and Tech. 43, 1 (2006) <https://www.tandfonline.com/doi/abs/10.1080/18811248.2006.9711062>`_

