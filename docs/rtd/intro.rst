************
Introduction
************

.. note::

   | Current Version: |version|
   | Date: |today|


Motivation and Physics Background
=================================

The fission of a heavy nucleus into two or more lighter fragments is usually accompanied with the emission of *prompt* neutrons and photons. In our current understanding of the fission process, the fission fragments are produced in a certain state of deformation and intrinsic excitation, eventually resulting in the production of excited fission fragments. Very quickly, those primary fragments will evaporate neutrons and photons to reach a more stable configuration, either a ground-state or a long-lived isomer. The post-neutron emission fission fragments, also called fission products, will possibly further :math:`\beta`-decay, leading to another burst of :math:`\beta`-delayed neutron and photon emissions.

The study of the prompt fission neutrons and photons is important to better model the nuclear fission process, constrain the collective and intrinsic configurations of the nascent fragments near the scission point, and understand the sharing of the available excitation energy between the two fragments. This study is also highly relevant for applications, ranging from nuclear energy safety and efficiency to non-proliferation and stockpile stewardship missions.

Until recently, most of the nuclear data evaluation work related to prompt fission neutrons and photons was limited to their average number, or multiplicity, and their average energy spectrum, as a function of incident neutron energy in the case of neutron-induced fission reactions. Even for those somewhat simple quantities, scarce experimental data exist only, limited to some important isotopes and incident neutron energies. Phenomenological models have been developed over the years, e.g., the so-called `Los Alamos model (LAM) <https://www.tandfonline.com/doi/abs/10.13182/NSE82-5>`_, mostly for calculating the average prompt fission neutron spectrum (PFNS). Evaluated data on prompt fission always originated from very scarce experimental data, leaving important gaps even in the most modern nuclear data libraries such as `ENDF/B-VIII.0 <https://www.nndc.bnl.gov/endf/b8.0/>`_.

The :program:`CGMF` code was developed to model the de-excitation of the fission fragments on an event-by-event basis, following the successive emissions of neutrons and photons. This is radically different from what had been done in the past, allowing an unprecedented level of predictions on distributions and correlations of neutrons, photons and fission fragments. The development of this code is accompanied by a host of new fission experiments that look at increasing levels of details and correlations among the vast quantity of fission data. Correlations and distributions of post-scission data will be very useful to constrain the still free parameters entering the :program:`CGMF` code and models.


Synopsis of the :program:`CGMF` code
====================================

The :program:`CGMF` code is based on two older codes developed at LANL: :program:`FFD` [LA-CC-10-003] and :program:`CGM` [LA-CC-11-018]. It performs Monte Carlo simulations of the decay of excited fission fragments by emission of prompt neutrons and gamma rays. The Hauser-Feshbach statistical theory of compound nuclear reactions is used to compute the emission probabilities at each step of the cascade. Monte Carlo histories are then recorded and analyzed.
 
The average prompt fission neutron multiplicity :math:`\overline{\nu}`, the prompt fission neutron multiplicity distribution :math:`P(\nu)`, the average prompt fission neutron multiplicity as a function of mass, charge and kinetic energy of the fragment, :math:`\overline{\nu}(A,Z,TKE)`, etc,  can all be extracted from :program:`CGMF` calculations. Similar quantities can also be obtained for prompt gamma rays. In addition, :math:`n-n`, :math:`n-\gamma`, and :math:`\gamma-\gamma` correlations can be studied both in energy and angle.


The :program:`FFD` code was the first fission fragment evaporation code to be developed, and used the Weisskopf-Ewing approximation to evaporate neutrons. Characteristics of the emitted prompt neutrons could be retrieved, but only little information could be inferred for the prompt :math:`\gamma`-ray data, and no specific discrete transitions in particular fragments could be tagged. Meanwhile, the :program:`CGM` code was being developed as a general Monte Carlo implementation of the traditional statistical nuclear reaction codes, e.g., GNASH, EMPIRE, TALYS, COH, using the very well-established Hauser-Feshbach statistical theory of nuclear reactions. With :program:`CGM`, one can follow the decay of an excited compound nucleus by evaporation of photons, neutrons, and light charged particles until it reaches its ground-state or a long-lived isomer. Monte Carlo histories could be followed one-by-one to study correlations and `exclusive` data.

Initially written in FORTRAN 95, significant parts of :program:`FFD` were re-written into C++ classes to work seamlessly with the C++ code :program:`CGM`, leading to the release of the :program:`CGMF` code that applies the physics of :program:`CGM` to the de-excitation of fission fragments.


.. _moreinfo:

For more information
=====================

This online user manual is intended to become the main reference for :program:`CGMF`. The main reference to cite :program:`CGMF` is its official documentation published as `Talou, Stetcu, Jaffke, Rising, Lovell and Kawano, submitted to Comp. Phys. Comm. (2020) <>`_. In addition, several :ref:`publications` and presentations might be of interest to the reader wanting more information on how the code is actually used for practical studies. 

