/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file smatrix.cpp

  \brief S-matrix elements

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "optical.h"

static inline double
lagrange(double h,double a,double b,double c,double e,double f, double g){
  return( ((g-a)/60+0.15*(b-f)+0.75*(e-c))/h );
}


/***********************************************************/
/*      S Matrix Element for Spherical Calculation         */
/***********************************************************/
Complex omSmatrix(int m, double mesh, double rm, Complex *d1, Complex *d2,
                  Wavefunc *wfn)
{
  Complex win,dwin,f,u,s;

  win.real = wfn->internal[m-3].real;
  win.imag = wfn->internal[m-3].imag;
  dwin.real = lagrange(mesh,wfn->internal[m-6].real,wfn->internal[m-5].real,
                            wfn->internal[m-4].real,wfn->internal[m-2].real,
                            wfn->internal[m-1].real,wfn->internal[m  ].real);
  dwin.imag = lagrange(mesh,wfn->internal[m-6].imag,wfn->internal[m-5].imag,
                            wfn->internal[m-4].imag,wfn->internal[m-2].imag,
                            wfn->internal[m-1].imag,wfn->internal[m  ].imag);
  f = rational(dwin.real,dwin.imag, win.real, win.imag);

  f.real *= rm;
  f.imag *= rm;

  f.real = f.real-d1->real;
  double si1 = f.imag+d1->imag;
  double si2 = f.imag-d1->imag;
  u = rational(f.real,si1,f.real,si2);
  s.real = u.real*d2->real-u.imag*d2->imag;
  s.imag = u.real*d2->imag+u.imag*d2->real;

  return(s);
}
