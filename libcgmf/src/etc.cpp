/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file etc.cpp

  \brief Miscellaneous functions

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>

using namespace std;

#include "physics.h"
#include "etc.h"
#include "maths.h"

double absolute(Complex *x)
{
  return( x->real * x->real + x->imag * x->imag );
}

Complex compsqrt(Complex *x)
{
  Complex c;

  double r = sqrt(x->real * x->real + x->imag * x->imag);
  double t = 0.5*acos(x->real/r);
  r=sqrt(r);

  c.real = r*cos(t);
  c.imag = r*sin(t);

  return(c);
}

Complex comppow(Complex *x, double p)
{
  Complex c;

  if(p==1.0) c = *x;
  else{
    double r = sqrt(x->real * x->real + x->imag * x->imag);
    double t = p*acos(x->real/r);
    r=pow(r,p);

    c.real = r*cos(t);
    c.imag = r*sin(t);
  }
  return(c);
}


Complex rational(double a1, double a2, double b1, double b2)
{
  Complex c;

  double d = b1*b1+b2*b2;
  c.real = (a1*b1+a2*b2)/d;
  c.imag = (a2*b1-a1*b2)/d;

  return(c);
}

#ifndef HAVE_MINMAX
int min(int x, int y){ return( (x<y) ? x : y); }
int max(int x, int y){ return( (x>y) ? x : y); }
#endif

double cfmin(double x, double y){ return( (x<y) ? x : y); }
double cfmax(double x, double y){ return( (x>y) ? x : y); }


double dep_to_jvol(double anum, double v, double r, double a)
{
  if(v==0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 4.0*PI/3.0*r*r*r;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2+1.0);

  return(v*x3);
}

double dep_to_jsrf(double anum, double w, double r, double a)
{
  if(w==0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 16.0*PI/anum*rx*rx*a;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2/3.0+1.0);

  return(w*x3);
}

double jvol_to_dep(double anum, double jv, double r, double a)
{
  if(jv==0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 4.0*PI/3.0*r*r*r;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2+1.0);

  return(jv/x3);
}

double jsrf_to_dep(double anum, double jw, double r, double a)
{
  if(jw==0.0) return(0.0);

  double rx = r*pow(anum,1.0/3.0);
  double x1 = 16.0*PI/anum*rx*rx*a;
  double x2 = PI*a/rx;
  double x3 = x1*(x2*x2/3.0+1.0);

  return(jw/x3);
}

double gaussian_weight(double ew, double ei, double ex)
{
  return( exp(-(ei-ex)*(ei-ex)/(2*ew*ew))/(ew*sqrt(PI2)) );
}


/**********************************************************/
/*     Laguerre Polynomial Function                       */
/**********************************************************/
double laguerre(int n, double a, double x)
{
  double p2=0.0;

  double p0 = 1.0;       if(n==0) return(p0);
  double p1 = 1.0+a-x;   if(n==1) return(p1);

  for(int i=2 ; i<=n ; i++){
    p2=((-x+2*i+a-1.0)*p1-(i+a-1.0)*p0)/(double)i;
    p0=p1; p1=p2;
  }

  return(p2);
}

/**********************************************************/
/*     Gamma Function                                     */
/**********************************************************/
double loggamma(double x)
{
  double v, w;

  v = 1.0;
  while (x < 8) {  v *= x;  x++;  }
  w = 1.0 / (x * x);
  return (((((((-0.0295506535947712  * w + 0.0064102564102564) * w
                -0.0019175269175269) * w + 0.0008417508417508) * w
                +0.0005952380952381) * w + 0.0007936507936508) * w
                -0.0027777777777778) * w + 0.0833333333333333) / x
                + 0.5 * LOG_2PI - log(v) - x + (x - 0.5) * log(x);
}

double gam(double x)
{
  if(x < 0) return( PI / (sin(PI * x) * exp(loggamma(1 - x))) );
            return( exp(loggamma(x)) );
}

/**********************************************************/
/*     Legendre Function                                  */
/**********************************************************/
double legendre(int n,double x)
{
  double p=0.0,p3=0.0;

  double p0 = 1.0;
  double p1 = x;
  double p2 = (3*x*x-1)/2;

  if(n==0)        p=p0;
  else if(n==1)   p=p1;
  else if(n==2)   p=p2;
  else{
    for(int i=2 ; i<n ; i++){
      double pn = (double)i;
      p3 = ((2*pn+1)*x*p2-pn*p1)/(pn+1);
      p1=p2;  p2=p3;
    }
    p=p3;
  }
  return(p);
}

double legendre1(int n,double x)
{
  double d = sqrt(1.0-x*x);
  if(n==0 || d==0.0) return(0.0);
  else    return( n*(-x*legendre(n,x)+legendre(n-1,x))/d );
}

/**********************************************************/
/*   Associated Legendre Polynomials   P^m_L              */
/**********************************************************/
double assocLegendrePol(int l, int m, double x)
{
  double p1 = 1.0;
  double p2 = sqrt((1.0-x)*(1.0+x));
  double p3 = 1.0;

  for(int i=1 ; i<=m ; i++){
    p1 *= -p3*p2;
    p3 += 2.0;
  }
  if(l==m) return(p1);

  p2 = x*(2*m+1)*p1;
  if(l==m+1) return(p2);

  for(int i=m+2 ; i<=l ; i++){
    p3 = (x*(2*i-1)*p2-(i+m-1)*p1)/(i-m);  p1 = p2;  p2 = p3;
  }
  return(p3);
}


/**********************************************************/
/*     Bessel Function I_{n+1/2}(x)                       */
/**********************************************************/
static const double EPSBESS = 1.0e-10;
double bessi2(int n, double x)
{
  double bi=0.0,x2,am,an,g1,g2,g3,s;
  int m,k;

  if(x<=1.0){
    x2=0.25*x*x;
    s=g1=1.0/(n+0.5);
    for(k=1;;k++){
      g1 *= x2/( k * (n+1+k*0.5) );
      s  += g1;
      if(fabs(g1)<EPSBESS) break;
    }
    bi=pow(0.5*x,(double)n+0.5)*s/gam(n+0.5);
  }
  else{
    if(     x<= 10.0) { am = 1.6  *x+11;  an = 1.1  *x+ 4; }
    else if(x<= 40.0) { am = 0.65 *x+18;  an = 0.466*x+11; }
    else if(x<=100.0) { am = 0.367*x+31;  an = 0.25 *x+21; }
    else              { am = 0.28 *x+39;  an = 0.2  *x+26; }
    m = (int)(am + ( ((double)n>an) ? n-an : 0.0 ));
    x2=2.0/x;
    s=g3=0.0; g2=1.e-36;
    for(k=m;k>=0;k--){
      g1 = (k+0.5)*x2*g2+g3;  if( k==n ) bi=g2;
      s += (k+0.5)*g2;
      g3=g2;
      g2=g1;
    }
    bi *= sqrt(x/(2*PI))*exp(x)/s;
  }

  return(bi);
}


/**********************************************************/
/*     Hankel Function K_{n+1/2}(x)                       */
/**********************************************************/
double bessk2(int n, double x)
{
  double bk,x2,g,s;
  int k;

  x2=1.0/(2*x);
  s=g=1.0;
  for(k=1;k<=n;k++){
    g *= x2*(n+k)*(n-k+1)/(double)k;
    s += g;
  }
  bk=sqrt(PI*x2)*s*exp(-x);
  return(bk);
}

