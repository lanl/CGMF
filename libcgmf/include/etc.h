/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file etc.h

  \brief Prototype of miscellaneous functions

*/


#ifndef __ETC_H__
#define __ETC_H__

double  jvol_to_dep           (double, double, double, double);
double  jsrf_to_dep           (double, double, double, double);
double  dep_to_jvol           (double, double, double, double);
double  dep_to_jsrf           (double, double, double, double);

#ifndef HAVE_MINMAX
int     min                   (int   ,int   );
int     max                   (int   ,int   );
#endif

double  cfmin                 (double,double);
double  cfmax                 (double,double);

double  gaussian_weight       (double, double, double);
double  laguerre              (int,    double, double);
double  gam                   (double);
double  loggamma              (double);
double  legendre              (int, double);
double  legendre1             (int, double);
double  assocLegendrePol      (int,int,double);
double  bessi2                (int, double);
double  bessk2                (int, double);

#endif //__ETC_H__
