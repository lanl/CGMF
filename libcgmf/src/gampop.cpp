/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file gampop.cpp

  \brief Gamma-ray population

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "cgm.h"
#include "mchf.h"
#include "global.h"
#include "config.h"

static double specTc2c (Statcalcmode, int, int, int, int, int, double, double, Nucleus *);
static double specTc2d (Statcalcmode, int, int, int, int, int, double, double, Nucleus *);


/**********************************************************/
/*      Gamma-ray Transition                              */
/**********************************************************/
double specTransitionGamma(Statcalcmode mode, int k0, int p0, int j0,
                           double **tg, double q, Nucleus *n, double *spc)
{
  /*** continuum to continuum */
  double x,sum = 0.0;
  for(int k1=k0+1 ; k1<n->ncont ; k1++){

    x   = specTc2c(mode,-1,k1,p0,j0,E1_MULTIPOL,tg[E1][k1],q,n);
    x  += specTc2c(mode, 1,k1,p0,j0,M1_MULTIPOL,tg[M1][k1],q,n);
    x  += specTc2c(mode, 1,k1,p0,j0,E2_MULTIPOL,tg[E2][k1],q,n);

    if(mode == hauser) spc[k1-k0] += x;
    sum += x;
  }

  /*** continuum to level */
  for(int i=0 ; i<n->ndisc ; i++){
    double eg = n->excitation[k0] - n->lev[i].energy;
    double ex = n->lev[i].energy;
    double a  = ldDensityParameter(ex,(double)n->za.getA(),&n->ldp);

    double te1 = gdrGammaTransmission(GL,E1_MULTIPOL,eg,E1,a  ,ex );
    double tm1 = gdrGammaTransmission(SL,M1_MULTIPOL,eg,M1,0.0,0.0);
    double te2 = gdrGammaTransmission(SL,E2_MULTIPOL,eg,E2,0.0,0.0);

    x   = specTc2d(mode,-1,i,p0,j0,E1_MULTIPOL,te1,q,n);
    x  += specTc2d(mode, 1,i,p0,j0,M1_MULTIPOL,tm1,q,n);
    x  += specTc2d(mode, 1,i,p0,j0,E2_MULTIPOL,te2,q,n);

    if(mode == hauser){
#ifdef HAVE_PRIVATE_ENERGY_GRID
      int kg = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,eg,n->binwidth);
#else
      int kg = specFindEnergyBin(eg,n->de);
#endif
      if(kg >= 0) spc[kg] += x;
    }

    sum += x;
  }

  return(sum);
}

 
/**********************************************************/
/*      Gamma-ray Transition : Continuum to Continuum     */
/**********************************************************/
double specTc2c(Statcalcmode mode, int pg, int k, int p0, int j0, int l,
                double tg, double q, Nucleus *n)
{
  double  x0=0.0,x1=0.0,dx=0.0,sum=0.0;
  int j1min = abs(j0-2*l);
  int j1max =     j0+2*l ;

  if(j1max > 2*MAX_J-1) j1max=2*MAX_J-1;
  if(n->jmax < j1max/2) n->jmax=j1max/2;

  for(int j1=j1min ; j1<=j1max ; j1+=2){
    if(j0==0 && j1==0) continue; // gamma decay prohibited for 0->0
    int idx = (j1-(int)(2.0*halfint(n->lev[0].spin)))/2;

    if(pg*p0==1){
      dx = n->density[k][idx].even*n->binwidth[k];
      x0 = tg*dx;
      if(mode==hauser) n->pop[k][idx].even += (x1=x0*q);
    }else{
      dx = n->density[k][idx].odd *n->binwidth[k];
      x0 = tg*dx;
      if(mode==hauser) n->pop[k][idx].odd  += (x1=x0*q);
    }
    if(mode==sumall) sum += x0;
    else if(mode == hauser) sum += x1;
  }

  return(sum);
}


/**********************************************************/
/*      Gamma-ray Transition : Continuum to Levels        */
/**********************************************************/
double specTc2d(Statcalcmode mode, int pg, int i, int p0, int j0, int l,
                double tg, double q, Nucleus *n)
{
  double x1=0.0,sum=0.0;
  double s1min = abs(j0/2.0-l);
  double s1max =     j0/2.0+l ;

  if(j0==0 && (int)n->lev[i].spin==0 && j0==0) return (sum);

  if(n->lev[i].spin >=s1min && n->lev[i].spin<=s1max){
    if(pg == p0*n->lev[i].parity){
      if(mode == hauser){
        n->lpop[i] += (x1=tg*q);
      }
      if(mode==sumall) sum += tg;
      else if(mode == hauser) sum += x1;
    }
  }

  return(sum);
}
