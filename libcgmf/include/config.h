/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file config.h

  \brief Customize CGMF, define data directory, system-dependent parameters

*/

#ifndef __CONFIG_H__
#define __CONFIG_H__

/*----------------------------------------------------------
Top directory in which CGM can find data files - 
     RIPL3 level data   /dir/to/levels/z001.dat
     kcksyst.dat        /dir/to/kcksyst.dat
*/

extern std::string datadir;


/*----------------------------------------------------------
Neutron emission channel
     This value allows / suppresses competing neutron emission
*/

const bool INCLUDE_NEUTRON_EMISSION = true;

// if excitation energy is Bn + this, gamma-ray emission ignored
// set zero to include all gammas

const double IGNORE_CAPTURE_GAMMA_RAY = 0.0;


/*----------------------------------------------------------
Energy bin width in the continuum
*/

//const double ENERGY_BIN = 0.10;  // 100 keV
//const double ENERGY_BIN = 0.01;  // 10 keV
//const double ENERGY_BIN = 0.05;  // 50 keV
extern double ENERGY_BIN;

/*----------------------------------------------------------
 Experimental time coincidence window (in seconds)
 By default, set to -99.0 to let all levels decay to the 
 ground.state.
 */

const double EXPERIMENTAL_TIME_WINDOW = -99.0; // [default]
// const double EXPERIMENTAL_TIME_WINDOW = 1.0e-8; // 10ns


/*----------------------------------------------------------
Beta strength function Gaussian broadening resolution
*/

const double PROFILE_BROADENING_WIDTH = 0.1;  // 100 keV


/*----------------------------------------------------------
Binary reaction calculation
     suppress multiple step reactions,
     calculate decay from the initial state only
*/

const bool BINARY_REACTION_ONLY = false;


/*----------------------------------------------------------
 Neutron spectrum calculated by Evaporation model
 no gamma-ray emission
 */
const bool EVAPORATION_SPECTRUM = false;


/*----------------------------------------------------------
Custom energy grid for spectra
*/

#undef HAVE_PRIVATE_ENERGY_GRID
//#define HAVE_PRIVATE_ENERGY_GRID

#ifdef   HAVE_PRIVATE_ENERGY_GRID
const int NUMBER_OF_PRIVATE_GRID = 272;
#define   GRID_STRUCTURE_FILE    "privategrid.h"
#endif


// no gamma decay if excitation energy of continuum is less than this
//const double CONTINUUM_LOWER_CUT = 0.02;
extern double CONTINUUM_LOWER_CUT;


/*----------------------------------------------------------
Discrete gamma-ray Internal Conversion control
     ICC considered in discrete transisions
*/

const bool INCLUDE_INTERNAL_CONVERSION = true;


/*----------------------------------------------------------
 Monte Carlo control
 RANDOM_SEED_BY_TIME
 change random number seed each time by using
 internal clock
 
 PERTURB_EXCITATION ENERGY
 random noise is added to the energies for
 the continuum to continuum transitions
 
 EVENT_OUTPUT_FORMAT
 0: output MC history in the CGM standard format
 1: output MC history in the GEANT4 input format
 2: output MC history in the FFD input format
 */

const bool RANDOM_SEED_BY_TIME = false;
const bool PERTURB_EXCITATON_ENERGY = true;
const int  EVENT_OUTPUT_FORMAT = 0;

#endif //__CONFIG_H__
