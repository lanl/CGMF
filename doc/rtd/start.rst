
===============
Getting Started
===============

Obtaining CGMF
--------------

CGMF is an open source code and can be cloned from the github repository at `https://github.com/lanl/CGMF <https://github.com/lanl/CGMF>`_ by typing:


.. code-block:: bash

	git clone git@github.com:lanl/CGMF.git


Installing CGMF
---------------

The build system for CGMF is based on the `CMake <https://cmake.org/>`_ build system generation tool. The current version requires CMake 3.16.2 minimumm and can be freely downloaded at `https://cmake.org/download/ <https://cmake.org/download>`_. Assuming that the :program:`CGMF` source tree is located in the $CGMFPATH directory, the Linux commands to build the code library and utilities are:

.. code-block:: console

	cd $CGMFPATH
	mkdir build
	cd build
	cmake ..
	make


This creates the static library `libcgmf.a` in the $CGMFPATH/build/libcgmf directory and also creates the executable `cgmf.x` in the $CGMFPATH/build/utils/cgmf directory.  Other options for building the code include:

.. code-block:: console

	* CMAKE\_BUILD\_TYPE=(Debug, RelWithDebInfo, Release [default])
	* CMAKE\_INSTALL\_PREFIX=(/usr/local/ [default])
	* cgmf.shared\_library=(ON, OFF [default])
	* cgmf.x.MPI=(ON, OFF [default])


Running CGMF
------------

To launch a CGMF run, type:

.. code-block:: c++

	./cgmf.x -i 98252 -e 0.0 -n 1000000

which means that CGMF is run for the spontaneous fission (incident energy is set to 0.0 (using `-e 0.0`) of 252Cf (ZAID=98252) with 1,000,000 fission events.


.. list-table:: CGMF Arguments
   :widths: 20 20 160
   :header-rows: 0

   * - -i $ZAIDt
     - [required]
     - 1000*Z+A of target nucleus, or fissioning nucleus if spontaneous fission

   * - -e $Einc
     - [required]
     - incident neutron energy in MeV (0.0 for spontaneous fission)

   * - -n $nevents
     - [required]
     -	number of Monte Carlo fission events to run or to be read.  If $nevents is negative, produces initial fission fragment yields Y(A,Z,KE,U,J,p)

   * - -t $timeCoinc
     - [optional]
     - time coincidence window for long-lived isomer gamma-ray emission cutoff (in sec)

   * - -d $datapath
     - [optional]
     -	overrides the environment variable CGMFDATA and default datapath

   * - -f $filename
     - [optional]
     - fission histories results file ("results.cgmf" is default)


The CGMF run above would create a history file (`histories.cgmf') as well as a concise summary of important average quantities on the console, such as:

.. code-block:: console

 	//// CGMF Results ////

	Reaction: spontaneous fission of (98,252)

	Average Light Fragment (Z,A) = (42.56,108.40)
	Average Heavy Fragment (Z,A) = (55.44,143.60)

	Average Kinetic Energies: LF = 105.73 MeV ; HF = 80.07 MeV ; <TKE> = 185.79 MeV
	Average Excitation Energies: LF = 18.33 MeV ; HF = 13.61 MeV ; <TXE> = 31.94 MeV

	Average Fragment Spins: <J>_LF = 9.11 hbar ; <J>_HF = 9.92 hbar ; <J> = 9.52 hbar

	*** Prompt Fission Neutrons ***

	Multiplicities (n/f):  <nu>_LF = 2.12 ; <nu>_HF = 1.69 ; <nu>_prefission = 0.00 ; <nu>_tot = 3.82 
	c-o-m Energies:  <Ecm>_LF = 1.34 MeV ; <Ecm>_HF = 1.21 MeV ; <Ecm>_prefission = 0.00 MeV ; <Ecm>_tot = 1.28 MeV
	Lab. Energies:   <Elab>_LF = 2.27 MeV ; <Elab>_HF = 1.72 MeV ; <Elab>_prefission = 0.00 MeV ; <Elab>_tot = 2.02 MeV

	*** Prompt Fission Gammas ***

	Multiplicities (g/f):  <nu_g>_LF = 4.30 ; <nu_g>_HF = 4.07 ; <nu_g>_tot = 8.37 
	Gamma Energies:   <Eg>_LF = 0.76 MeV ; <Eg>_HF = 0.74 MeV ; <Eg>_tot = 0.75 MeV

	//// THE END ////



