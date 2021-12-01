/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file exciparm.cpp

  \brief Parameter calculation used in exciton model

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <algorithm>

using namespace std;

#include "exciton.h"

/*** log(N!) */
const double fact[] = 
 {0.00000e+00,  0.00000e+00,  6.93147e-01,  1.79176e+00,  3.17805e+00,
  4.78749e+00,  6.57925e+00,  8.52516e+00,  1.06046e+01,  1.28018e+01,
  1.51044e+01,  1.75023e+01,  1.99872e+01,  2.25522e+01,  2.51912e+01,
  2.78993e+01,  3.06719e+01,  3.35051e+01,  3.63954e+01,  3.93399e+01,
  4.23356e+01,  4.53801e+01,  4.84712e+01,  5.16067e+01,  5.47847e+01,
  5.80036e+01,  6.12617e+01,  6.45575e+01,  6.78897e+01,  7.12570e+01,
  7.46582e+01};


/**********************************************************/
/*      Single Particle State Density                     */
/**********************************************************/
double preqSingleStateDensity(int x)
{
  /*

   WARNING: This parameter is only for neutron-induced reactions
   */
  const double f = 19.2;
  double g = x/f;
  
  return(g);
}


/**********************************************************/
/*      Pairing Energy - Delta                            */
/**********************************************************/
double preqPairingDelta(int z, int a)
{
  int n = a-z;
  int c = 1;
  if(      (n%2)==1 && (z%2)==1 ) c = 0;
  else if( (n%2)==0 && (z%2)==0 ) c = 2;
  else c = 1;
  return( (double)c*12.0 / sqrt((double)a) );
}


/**********************************************************/
/*      Fu's Pairing Correction Factor                    */
/**********************************************************/
double preqPairingCorrection(int n, double ex, SPdens *p)
{
  if(p->pairing == 0.0) return(0.0);

  double x  = 0.0;
  double r  = ex/p->pairing;
  double g  = p->gz+p->gn;
  double nc = 1.38629436111989 * g * sqrt(p->pairing/g)*8.0/7.0;

  if( r >= (0.716 + 2.44*pow((double)n/nc,2.17))){
    x = 0.996 - 1.76*pow((double)n/nc,1.6) * pow(r,-0.68);
    x = p->pairing*(1.0-x*x);
  }
  else{
    x = p->pairing;
  }

  return(x);
}


/**********************************************************/
/*      Pauli Blocking Correction                         */
/**********************************************************/

double preqPauliCorrection(int zp, int zh, int np, int nh, double gz, double gn)
{
  if(zp<0 || zh<0 || np<0 || nh<0) return(0.0);

  double pm = (double)max(zp,zh);
  double nm = (double)max(np,nh);
  double cp = ( pm*pm - (zp*zp + zh*zh + zp+zh)/4.0 )/gz;
  double cn = ( nm*nm - (np*np + nh*nh + np+nh)/4.0 )/gn;
  return(cp+cn);
}


/**********************************************************/
/*      Averaged Matrix Elements by Koning                */
/**********************************************************/
double preqMSquare(int a, int n, Preeq *q)
{
  double a3   = pow(a,3.0);
  double msq  = ( 6.8 + 4.2e+5/pow(q->ex_total/(double)n+10.7,3.0) )/a3;
  double msqx = 1.5 * msq;

  q->m2zz = q->m2nn = msq;
  q->m2zn = q->m2nz = msqx;

  return(msq);
}


/**********************************************************/
/*      Finite Potential Well Function                    */
/*      This includes Kalbach surface effect for h=1      */
/**********************************************************/
double preqFiniteWell(int n, int h, double ex, double v)
{
  const double ef = 38.0;
  if(h>1) v = ef;

  double x = 1.0;

  if(v<(ef-0.5)){
    double w = v*(ef-v)/(2.0*ef);

    if(n==1 && h==1){
      if(ex > 1.16*ef)  x = 0.0;
      else if(ex < v)   x = 1.0;
      else              x = 1.0/(1.0+exp( (ex-v)/w ));
    }
    else{
      double qsum = 0.0;
      double xsum = 0.0;
      for(int p = -4 ; p <=4 ; p++){
        double v1 = v + p*w;
        if(v1<0.0 || v1>ef) continue;

        double q  = 1.0/(1.0+exp( (v-v1)/w ))/(1.0+exp( (v1-v)/w ));

        double x1 = 1.0;
        for(int i=1 ;  i<=h ;i++){
          if(ex - i*v1>0) x1 += ((i%2==0) ? 1: -1)
                               * exp(fact[h]-fact[i]-fact[h-i])
                               * pow(1.0-i*v1/ex,(double)n-1.0);
        }
        qsum += q;
        xsum += q*x1;
      }
      if(qsum>0.0) x = xsum / qsum;
      if(x<0.0) x=0.0;
    }
  }else{
    if(ex>v){
      for(int i=1 ;  i<=h ;i++){
        if(ex - i*v>0) x += ((i%2==0) ? 1: -1)
                           * exp(fact[h]-fact[i]-fact[h-i])
                           * pow(1.0-i*v/ex,(double)n-1.0);
      }
    }
  }

  return(x);
}


/**********************************************************/
/*      Finite Potential Well Depth Parameter             */
/**********************************************************/
double preqPotentialDepth(double e, int zi, int at)
{
  double a3 = pow((double)at,-1.0/3.0);
  double e4 = pow(e,4.0);
  double v  = (zi==0) ?
    12.0+26.0*e4/(e4+pow(245.0*a3,4.0)) :
    22.0+16.0*e4/(e4+pow(450.0*a3,4.0)) ;
  return(v);
}


/**********************************************************/
/*      State Density by Betak and Dobes                  */
/**********************************************************/
double preqStateDensity(double ex, SPdens *spd, Exconf *exc)
{
  double omega = 0.0;

  int p  = exc->zp + exc->np;
  int h  = exc->zh + exc->nh;
  int n  = p  + h ;

  if(n==0) return(omega);
  if(exc->zp<0 || exc->zh<0 || exc->np<0 || exc->nh<0) return(omega);
  if(ex<=0.0) return(omega);

  double ux = ex - preqPairingCorrection(n,ex,spd);
  double ap = preqPauliCorrection(exc->zp,exc->zh,exc->np,exc->nh,spd->gz,spd->gn);

  if(ux<0.0 || ux<ap) return(omega);

  omega =  pow(spd->gz,exc->zp+exc->zh) * pow(spd->gn,exc->np+exc->nh)
          /(exp(fact[exc->zp]+fact[exc->zh]+fact[exc->np]+fact[exc->nh]+fact[n-1]))
          * pow((ux-ap),(double)n-1.0) 
          * preqFiniteWell(n,h,ux,spd->well_depth);

  return(omega);
}

