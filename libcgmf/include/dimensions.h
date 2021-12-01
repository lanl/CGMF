/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file dimensions.h

  \brief Array sizes & default parameters

*/

#ifndef __DIMENSIONS_H__
#define __DIMENSIONS_H__

/****************************/
/*       ARRAY SIZES        */
/****************************/

// array sizes defined in cgmfEvents.cpp at initialization of cgmfEvent or cgmEvent
extern int MAX_ENERGY_BIN;
extern int MAX_J;
extern int MAX_LEVELS;

const int MAX_MULTIPOL     =    3 ;  /* multipolarity, E1, E1(def), M1, E2    */
const int MAX_GDR          =    6 ;  /* maximum GDRs and pygmy resonance      */
const int MAX_CHANNEL      =    2 ;  /* maximum decay channel                 */
const int MAX_GAMMA_LINES  = 1000 ;  /* maximum line gamma-rays               */
const int MAX_COMPOUND     =   10 ;  /* maximum number of compound nucleus    */
const int SPECTRA_OUTPUT   =    4 ;  /* number of spectra to be calculated    */


/****************************/
/*     DEFAULT PARAMETERS   */
/****************************/

const double NORM_FACT      = 10.0 ; /* conversion factor to [mb]             */
const double SPIN_CUTOFF    = 1e-8 ; /* cut-off of total anuglar momentum     */


#endif //__DIMENSIONS_H__
