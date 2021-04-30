/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file intwave.cpp

  \brief Wave functions inside nucleus, Fox-Goodwin or Numerov method

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include "optical.h"
#include "etc.h"


/**********************************************************/
/*      Solve Schroedinger Equation for Optical Potential */
/**********************************************************/
void omInternalFunction(int integ, double delta, double wavesq,
                        double xl, double xs, double xj,
                        Potential *pot, Wavefunc *wfn)
{
  double xl_term,so_term,deltasq;
  Complex a1,a2,a3,q1,q2,q3,p1,p2,p3;
  Complex complex_Zero(0.0,0.0);

  deltasq = delta*delta/12.0;
  xl_term = xl*(xl+1.0);
  so_term = xj*(xj+1.0)-xl_term-xs*(xs+1.0);

/***************************************/
/*     Fox-Goodwin Method              */
/***************************************/
  wfn->internal[0]=p1=complex_Zero;
  wfn->internal[1].real=p2.real=pow(delta,xl+1.0);
  wfn->internal[1].imag=p2.imag=p2.real*(-1.0e-10);

  q1     =complex_Zero;
  q2.real=xl_term/pot->radi[1]-wavesq
          -pot->mean_field[1].real -so_term*pot->spin_orbit[1].real;
  q2.imag=-pot->mean_field[1].imag -so_term*pot->spin_orbit[1].imag;
   
  for(int i=2 ; i<=integ ; i++){
    q3.real=xl_term/pot->radi[i]-wavesq
            -pot->mean_field[i].real-so_term*pot->spin_orbit[i].real;
    q3.imag=-pot->mean_field[i].imag-so_term*pot->spin_orbit[i].imag;

    a1.real=1.0     -deltasq*q1.real;
    a1.imag=        -deltasq*q1.imag;
    a2.real=2.0+10.0*deltasq*q2.real;
    a2.imag=    10.0*deltasq*q2.imag;

    p3.real=a2.real*p2.real - a2.imag*p2.imag 
           -a1.real*p1.real + a1.imag*p1.imag;
    p3.imag=a2.real*p2.imag + a2.imag*p2.real 
           -a1.real*p1.imag - a1.imag*p1.real;

    a3=rational(1.0, 0.0, 1.0-deltasq*q3.real, -deltasq*q3.imag);

    wfn->internal[i].real=a3.real*p3.real-a3.imag*p3.imag;
    wfn->internal[i].imag=a3.real*p3.imag+a3.imag*p3.real;

    q1=q2; q2=q3;
    p1=p2; p2=wfn->internal[i];
  }
}



/**********************************************************/
/*      Normalization Factor of  Wave Function            */
/**********************************************************/
Complex omNormalizationFactor(int l,int i,double yeta,double sig0,Complex *s,Wavefunc *wfn)
{
  Complex c,x,w,coul;
  double d;

  c.real=     s->imag /2.0;
  c.imag=(1.0-s->real)/2.0;

/* W = F + C(G+iF) */
  w.real=  wfn->external[l].real*c.real + wfn->external[l].imag*(1-c.imag);
  w.imag=  wfn->external[l].imag*c.real + wfn->external[l].real*   c.imag ;

  d     =  absolute(&(wfn->internal[i]));
  x.real= (w.real*wfn->internal[i].real+w.imag*wfn->internal[i].imag)/d;
  x.imag= (w.imag*wfn->internal[i].real-w.real*wfn->internal[i].imag)/d;

  if(yeta!=0.0){
    coul = omCoulombPhaseFactor(l,yeta,sig0);
    w = x;
    x.real = w.real*coul.real - w.imag*coul.imag;
    x.imag = w.real*coul.imag + w.imag*coul.real;
  }

  return(x);
}


/**********************************************************/
/*      Coulomb Phase Factor                              */
/**********************************************************/
Complex omCoulombPhaseFactor(int l,double yeta,double sig)
{
  Complex coul;

  for(int l0=1 ; l0<=l ; l0++) sig += atan(yeta/l0);
  coul.real=cos(sig);
  coul.imag=sin(sig);

  return(coul);
}


/**********************************************************/
/*      Normalize Wave Function (Distorted Wave)          */
/**********************************************************/
void omNormalization(int integ, double beta, Complex *norm, Potential *pot, Wavefunc *wfn)
{
  double alpha=beta*beta/4.0;
  double c=1.0;
  for(int i=0 ; i<=integ ; i++){
    Complex w=wfn->internal[i];
    if(beta>0.0){
      c=1.0/sqrt(1.0 + alpha*(pot->mean_field[i].real+pot->coulomb[i]));
    }
    wfn->internal[i].real=(norm->real*w.real-norm->imag*w.imag)*c;
    wfn->internal[i].imag=(norm->real*w.imag+norm->imag*w.real)*c;
  }
}

