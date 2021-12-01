/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file parpop.cpp

  \brief Particle population

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "cgm.h"
#include "config.h"

/*** transmission coefficient index calculator */
static inline int tj_index(int a, int b, int c){
  return((a*3 + (1-b)*c)/2);
}

static double  specTransitionParticleCont(Statcalcmode, int, int, int,
                                          Transmission *, double, Nucleus *, Nucleus *);
static double  specTransitionParticleDisc(Statcalcmode, int, int, int,
                                          Transmission *, double, Nucleus *, Nucleus *);


/**********************************************************/
/*      Particle Transition                               */
/**********************************************************/
double specTransitionParticle(Statcalcmode mode, int k0, int pcn, int jcn,
                              Transmission  *tc, Transmission  *td, double q,
                              Nucleus *n0, Nucleus *n1, double *spc)
{
  double x,sum = 0.0;
  Transmission *tp;

  /*** Continuum to continuum : loop over daughter excitation */
  for(int k1=k0 ; k1<n1->ncont ; k1++){
    int kp = k1-k0;
    tp = &tc[kp];
    if(tp->lmax<=0) continue;

    x = specTransitionParticleCont(mode,k1,pcn,jcn,tp,q,n0,n1);

    if(mode == hauser) spc[kp] += x;

    sum += x;
  }

  /*** Continuum to levels: Loop over daughter excitation */
  for(int k1=0 ; k1<n1->ndisc ; k1++){
    tp = &td[k1];
    if(tp->lmax<=0) continue;

    x = specTransitionParticleDisc(mode,k1,pcn,jcn,tp,q,n0,n1);

    if(mode == hauser){
      double ep = n0->excitation[k0]-n1->lev[k1].energy-n0->cdt[neutron].binding_energy;

#ifdef HAVE_PRIVATE_ENERGY_GRID
      int kp = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,ep,n1->binwidth);
#else
      int kp = specFindEnergyBin(ep,n1->de);
#endif
      if(kp >= 0) spc[kp] += x;
    }

    sum += x;
  }

  return(sum);
}


/**********************************************************/
/*      Particle Transition : Continuum to Continuum      */
/**********************************************************/
double specTransitionParticleCont(Statcalcmode mode, int k1, int pcn, int jcn,
                                  Transmission *tc, double q, Nucleus *n0, Nucleus *n1)
{
  int    jtop = 0;
  double x0=0.0,x1=0.0,dx=0.0,sum=0.0;

  /*** Loop over particle angular momentum */
  for(int lp=0 ; lp<=2*tc->lmax ; lp+=2){
    int pl = ((lp/2)%2==1) ? -1 : 1;
    for(int sp=n0->cdt[neutron].spin2 ; sp>=-n0->cdt[neutron].spin2 ; sp-=2){
      int jp = lp + sp;  if(jp<0) continue;

      double tj  = tc->tran[tj_index(lp,sp,n0->cdt[neutron].spin2)];

      int j1min = abs(jcn-jp);
      int j1max =     jcn+jp ; if(j1max>2*(MAX_J-1)) j1max=2*(MAX_J-1);

      if(j1max > jtop) jtop = j1max;

      for(int j1=j1min ; j1<=j1max ; j1+=2){

        int idx = (j1-(int)(2.0*halfint(n1->lev[0].spin)))/2;

        if(pl*pcn==1){
          dx = n1->density[k1][idx].even*n1->binwidth[k1];
          x0 = tj*dx;
          if(mode == hauser) n1->pop[k1][idx].even += (x1=q*x0);
        }else{
          dx = n1->density[k1][idx].odd*n1->binwidth[k1];
          x0 = tj*dx;
          if(mode == hauser) n1->pop[k1][idx].odd  += (x1=q*x0);
        }
        if(mode==sumall) sum += x0;
        else if(mode == hauser) sum += x1;
      }
    }
  }
  if(n1->jmax < jtop/2) n1->jmax = jtop/2;
  return(sum);
}


/**********************************************************/
/*      Particle Transition : Continuum to Levels         */
/**********************************************************/
double specTransitionParticleDisc(Statcalcmode mode, int k1, int pcn, int jcn,
                                  Transmission *td, double q, Nucleus *n0, Nucleus *n1)
{
  double x1=0.0,sum=0.0;

  int j1 = (int)(2.0*n1->lev[k1].spin);

  /*** Loop over particle angular momentum */
  for(int lp=0 ; lp<=2*td->lmax ; lp+=2){
    int pl = ((lp/2)%2==1) ? -1 : 1;
    if(pl != pcn*n1->lev[k1].parity) continue;

    for(int sp=n0->cdt[neutron].spin2 ; sp>=-n0->cdt[neutron].spin2 ; sp-=2){
      int jp = lp + sp;  if(jp<0) continue;

      int j1min = abs(jcn-jp);
      int j1max =     jcn+jp ;

      if(j1>=j1min && j1<=j1max){

        double tj  = td->tran[tj_index(lp,sp,n0->cdt[neutron].spin2)];

        if(mode == hauser) n1->lpop[k1] += (x1=tj*q);
        if(mode==sumall) sum += tj;
        else if(mode == hauser) sum += x1;
      }
    }
  }

  return(sum);
}


