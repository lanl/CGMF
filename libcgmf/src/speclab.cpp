/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file speclab.cpp

  \brief Convert CM spectrum into LAB

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "cgm.h"
#include "global.h"

#include "maths.h"

static const double eps = 1.0e-99;

static double cmsSpecInterpolate(double, double *, double *);
static inline double lowfilter (double x){if(fabs(x)<eps) return 0.0; else return(x);}

/**********************************************************/
/*      Main Spectra Calculation                          */
/**********************************************************/
void cgmLabSpectrum(double ef, double *de, double **spc)
{
  double *spl;

  ef  = sqrt(ef);
  spl = new double [MAX_ENERGY_BIN];

  double emin = 0.0;
  double emax = de[0]*0.5;
  double etop = 0.0;

  int k0 = cgmZeroCut(2,spc);
  for(int k=0 ; k<=k0 ; k++) etop += de[k];

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){

    double e  = sqrt((emin+emax)/2.0);
    double e1 = e - ef;   e1 = e1*e1;  // {sqrt(E) - sqrt(Ef)}^2
    double e2 = e + ef;   e2 = e2*e2;  // {sqrt(E) + sqrt(Ef)}^2
    double e3 = (e2+e1)*0.5;
    double e4 = (e2-e1)*0.5;

    if(e1 > etop){
      k0 = k;
      break;
    }

    /*** Gauss-Legendre Integration in [e1,e2] */
    spl[k] = 0.0;
    for(int ig=0 ; ig<MAX_GAUSS ; ig++){
      double x1 = e3 + e4*gauss_x[ig];
      double x2 = e3 - e4*gauss_x[ig];
      double s1 = cmsSpecInterpolate(x1,de,spc[1]);
      double s2 = cmsSpecInterpolate(x2,de,spc[1]);
      spl[k] += (s1/sqrt(x1) + s2/sqrt(x2)) * gauss_a[ig];
    }
    spl[k] *= e4*0.25/ef;  // int phi(e)/(4*sqrt(e Ef)) de

    emin = emax;
    emax += de[k+1];
  }

  cout << "#      Emin        Emax        Ecal      "<< endl;
  cout << "#                                        " << "  Neutron   " << endl;

  emin = 0.0;
  emax = de[0]*0.5;

  for(int k=0 ; k<=k0 ; k++){
    double e    = (emin+emax)/2.0;

    cout << setw(5) << k;
    cout << setprecision(4) << setiosflags(ios::scientific);
    cout << setw(12) << emin
         << setw(12) << emax
         << setw(12) << e;
    cout << setw(12) << lowfilter(spl[k]) << endl;

    emin = emax;
    emax += de[k+1];
  }

  delete [] spl;
  return;
}


/***********************************************************/
/*      Interpolate CMS Neutron Spectrum                   */
/***********************************************************/
double cmsSpecInterpolate(double x, double *de, double *s)
{
  double y = 0.0;

  if(x <= 0.0)  return(y);

  double e1min = 0.0;
  double e1max = de[0]*0.5;
  double e1    = (e1min + e1max)*0.5;

  double e2min = e1max;
  double e2max = e2min + de[1];
	double e2;


  /*** if lower than the first point, interpolate by sqrt(E) */ 
  if( x < e1 ){
    y = s[0] * sqrt(x/e1);
    return(y);
  }

  /*** linear interpolation */
  for(int k=1 ; k<MAX_ENERGY_BIN ; k++){
    e1    = (e1min + e1max)*0.5;
    e2    = (e2min + e2max)*0.5;

    if( (e1 <=x ) && (x < e2) ){
      y = (x-e1)/(e2-e1) * (s[k]-s[k-1]) + s[k-1];
      break;
    }

    e1min = e1max;
    e1max = e2max;
    e2min = e2max;
    e2max = e2min + de[k+1];
  }

  return(y);
}
