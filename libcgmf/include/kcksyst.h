/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file kcksyst.h

  \brief Koura-Chiba-Kawano systematics for Gilbert-Cameron-Ignatyuk level density parameters

*/
#ifndef __KCKSYST_H__
#define __KCKSYST_H__

#include "structur.h"

#define KCKDATAFILE "kcksyst.dat"

/**************************************/
/*      kcksyst.cpp                   */
/**************************************/
int     kckDataRead                     (ZAnumber *, LevelDensity *);
double  kckAsymptoticLevelDensity       (double);
double  kckSpinCutoff                   (double);
double  kckTemperature                  (double, double);
double  kckE0                           (double, double, double);

int readkcksystdat();
int getkcksystdat(ZAnumber *, LevelDensity *);

#endif //__KCKSYST_H__
