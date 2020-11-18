/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file omsetform.cpp

  \brief Calculate radial part of optical potentials

*/

#include <iostream>
#include <cmath>

using namespace std;

#include "physics.h"
#include "optical.h"
#include "terminate.h"

static Complex omPotentialRadialMean    (unsigned int, double, Optical *);
static Complex omPotentialRadialSpo     (double, Optical *);
static double  omWoodsSaxon             (double, double, double);
static double  omDerivWoodsSaxon        (double, double, double);
static double  omGaussian               (double, double, double);
static double  omThomas                 (double, double, double);


/**********************************************************/
/*     Spherical Optical Potential Radial Form            */
/**********************************************************/
int omPotentialForm(unsigned int potform, int zzprod, Optical *omp, CCdata *cdt, Potential *pot)
{
  double r,x[6];
  int    n;

  double eps=log(CRIT_MATCHING);
  x[0]=x[1]=x[2]=x[3]=x[4]=x[5]=0.0;

  if(potform & 0x0010) x[0]=omp->R0+omp->a0   *(log(    omp->volume.real  /cdt->energy)-eps);
  if(potform & 0x0020) x[1]=omp->R0s+omp->a0s *(log(4.0*fabs(omp->surface.real)/cdt->energy)-eps);
  if(potform & 0x0001) x[2]=omp->Rv+omp->av   *(log(    omp->volume.imag  /cdt->energy)-eps);
  if(potform & 0x0002) x[3]=omp->Rs+omp->as   *(log(4.0*omp->surface.imag /cdt->energy)-eps);
  if(potform & 0x0004) x[3]=omp->Rs+omp->as   *sqrt(log(4.0*omp->surface.imag /cdt->energy)-eps);

  if(omp->spin_orbit.real!=0.0) x[4]=omp->Rvso+omp->avso*(log(     omp->spin_orbit.real *CSPO
                                                        /(omp->avso*omp->Rvso*cdt->energy))-eps);
  if(omp->spin_orbit.imag!=0.0) x[5]=omp->Rwso+omp->awso*(log(fabs(omp->spin_orbit.imag)*CSPO
                                                        /(omp->awso*omp->Rwso*cdt->energy))-eps);

  pot->rad_match=x[0];
  for(int i=1 ; i<6 ; i++) if(x[i] > pot->rad_match) pot->rad_match=x[i];

  n=(int)(pot->rad_match/pot->width)+1;
  r=pot->width*n; if(r<pot->rad_match) n++;
  pot->rad_match=pot->width*n;

  n+=3;  pot->n_match = n;

  if(n>MAX_POINTS) cgmTerminateCode("integration number too large",n);

  omPotentialFixedLength(n,potform,zzprod,omp,cdt,pot);

  return(n);
}


/**********************************************************/
/*     Fadial Form for Fixed Matching Radius              */
/**********************************************************/
int omPotentialFixedLength(int n, unsigned int potform, int zzprod, Optical *omp, CCdata *cdt, Potential *pot)
{
  double r;

  pot->rho_match = pot->rad_match*sqrt(cdt->wavesq);

  double c1 = 2.0*cdt->reduced_mass*AMUNIT/VLIGHTSQ/HBARSQ;
  double c2 = c1*CSPO;
  double c3 = c1*zzprod*PERMITTIV*COULOMBSQ;

  for(int i=0 ; i<=n ; i++){
    r=(double)i*pot->width;
    pot->radi[i]=r*r;
    pot->mean_field[i] = omPotentialRadialMean(potform,r,omp);
    pot->spin_orbit[i] = omPotentialRadialSpo(r,omp);
    pot->coulomb[i]    = (zzprod!=0) ? omPotentialRadialCoulomb(r,omp) : 0.0;

    pot->mean_field[i].real *= c1;
    pot->mean_field[i].imag *= c1;

    pot->spin_orbit[i].real *= c2;
    pot->spin_orbit[i].imag *= c2;

    pot->coulomb[i]         *= c3;

    pot->mean_field[i].real -= pot->coulomb[i];
  }

  return(n);
}



double omPotentialRadialCoulomb(double r, Optical *omp)
{
  return( (r <= omp->Rc) ? (3.0-r*r/(omp->Rc*omp->Rc))/(2.0*omp->Rc) : 1.0/r );
}



Complex omPotentialRadialMean(unsigned int potform, double r, Optical *omp)
{
  Complex p(0.0,0.0);
  unsigned int real_pot,imag_pot;

  potform = potform & 0x00ff;
  real_pot= potform >>     4;
  imag_pot= potform & 0x000f;

  switch(real_pot){
  case 1 : p.real =omp->volume.real *      omWoodsSaxon(r,omp->R0 ,omp->a0 );  break;
  case 2 : p.real =omp->surface.real* omDerivWoodsSaxon(r,omp->R0s,omp->a0s);  break;
  case 3 : p.real =omp->volume.real *      omWoodsSaxon(r,omp->R0 ,omp->a0 )
                  +omp->surface.real* omDerivWoodsSaxon(r,omp->R0s,omp->a0s);  break;
  default: p.real=0.0;
  }

  switch(imag_pot){
  case 1 : p.imag =omp->volume.imag       *omWoodsSaxon(r,omp->Rv,omp->av);    break;
  case 2 : p.imag =omp->surface.imag* omDerivWoodsSaxon(r,omp->Rs,omp->as);    break;
  case 3 : p.imag =omp->volume.imag       *omWoodsSaxon(r,omp->Rv,omp->av)
                  +omp->surface.imag* omDerivWoodsSaxon(r,omp->Rs,omp->as);    break;
  case 4 : p.imag =omp->surface.imag        *omGaussian(r,omp->Rs,omp->as);    break;
  case 5 : p.imag =omp->volume.imag       *omWoodsSaxon(r,omp->Rv,omp->av)
                  +omp->surface.imag        *omGaussian(r,omp->Rs,omp->as);    break;
  default: p.imag=0.0;
  }

  return(p);
}


Complex omPotentialRadialSpo(double r, Optical *omp)
{
  Complex p(0.0,0.0);
  if(omp->spin_orbit.real!=0.0) p.real =omp->spin_orbit.real*omThomas(r,omp->Rvso,omp->avso);
  if(omp->spin_orbit.imag!=0.0) p.imag =omp->spin_orbit.imag*omThomas(r,omp->Rwso,omp->awso);
  return(p);
}



double omWoodsSaxon(double x, double r, double a)
{
  return(  1/( 1+exp((x-r)/a) )  );
}

double omDerivWoodsSaxon(double x, double r, double a)
{
  double y = exp((x-r)/a);
  return(  4*y/( (1+y)*(1+y) )  );
}

double omGaussian(double x, double r, double a)
{
  double y = (x-r)*(x-r)/(a*a);
  return(  exp(-y)  );
}

double omThomas(double x, double r, double a)
{
  if(x==0.0) return(0.0);
  double  y=exp((x-r)/a);
  return(  y/( x*a*(1+y)*(1+y) )  );
}

