=============
The CGMF Code
=============

:program:`CGMF` is a code that simulates the emission of prompt fission neutrons and gamma rays from excited fission fragments right after scission. It implements a Monte Carlo version of the Hauser-Feshbach statistical theory of nuclear reactions to follow the decay of the fission fragments on an event-by-event basis. Probabilities for emitting neutrons and gamma rays are computed at each stage of the decay. Each fission event history records characteristics of the parent fragment (mass, charge, kinetic energy, momentum vector, excitation energy, spin, parity) and the number (multiplicity) and characteristics (energy, direction) of the prompt neutrons and gamma rays emitted in this event.

.. admonition:: Recommended publication for citing
   :class: tip

   Patrick Talou, Ionel Stetcu, Patrick Jaffke, Michael E. Rising, Amy E. Lovell, and Toshihiko Kawano, "Fission Fragment Decay Simulations with the CGMF Code," to be submitted to Comp. Phys. Comm. (2020).

.. admonition:: Support

   For any questions related to :program:`CGMF`, its use, and its code source, please email us at: :email:`cgmf-help@lanl.gov <cgmf-help@lanl.gov>`.

.. only:: html

--------
Contents
--------

.. toctree::
   :maxdepth: 2

   Introduction <intro>
   Getting Started<start>
   Physics Models<models>
   Code Details<code>
   Example Jupyter Notebooks<examples_notebook>
   Publications<publications>
   License Agreement<license>

.. note::

	| Current Code Version: |version|
	| Date: |today|

