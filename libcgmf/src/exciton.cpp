/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file exciton.cpp

  \brief Exciton model for precompound reactions

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "physics.h"
#include "structur.h"
#include "excinterface.h"
#include "exciton.h"

static double  preqEmissionRate     (Nucleus *, Pdata *,Transmission **, Preeq *,Exconf *, double **);
static void    preqLifeTime         (int, int, double **, double **, double **,
                                     double **, double **, double **, double **,
                                     Preeq *, Exconf *);
static inline double preqTransmissionAverage(int, int, double *);


/**********************************************************/
/*      Exciton Model for Pre-Equilibrium Reaction        */
/**********************************************************/
void    preqExcitonModel(Nucleus *ncl, System *sys, Pdata *pdt, Transmission **tc, double **spc)
{
  Preeq    prq;               // preequilibrium parameters
  Exconf   exc;               // exciton configuration
  SPdens   spd[MAX_CHANNEL];  // single-particle state densities

  double   *wkc[MAX_CHANNEL]; // continuum particle emission rate
  double   *pop[MAX_PREEQ];   // occupation probability
  double   *wk[MAX_PREEQ];    // integrated emission rate
  double   *tau[MAX_PREEQ];   // lifetime of exciton state
  double   *tap[MAX_PREEQ];   // lifetime only for exchange interaction
  double   *tzp[MAX_PREEQ];   // lambda z+
  double   *tnp[MAX_PREEQ];   // lambda n+
  double   *tz0[MAX_PREEQ];   // lambda z0
  double   *tn0[MAX_PREEQ];   // lambda n0

//---------------------------------------
//      Memory Allocation
  for(int j=0 ; j<MAX_CHANNEL ; j++){
    wkc[j] = new double [MAX_ENERGY_BIN];
  }
  for(int j=0 ; j<MAX_PREEQ   ; j++){
    pop[j] = new double [MAX_PREEQ];
    tzp[j] = new double [MAX_PREEQ];
    tnp[j] = new double [MAX_PREEQ];
    tz0[j] = new double [MAX_PREEQ];
    tn0[j] = new double [MAX_PREEQ];
    wk [j] = new double [MAX_PREEQ];
    tau[j] = new double [MAX_PREEQ];
    tap[j] = new double [MAX_PREEQ];
  }

  /*** Indices of f(zp,zh,np,nh) are hole numbers, 
       assuming that initial hole number is zero */
  for(int i=0 ; i<MAX_PREEQ ; i++){
    for(int j=0 ; j<MAX_PREEQ ; j++){
      pop[i][j] = wk[i][j]  = tau[i][j] = tap[i][j] = 0.0;
      tzp[i][j] = tnp[i][j] = tz0[i][j] = tn0[i][j] = 0.0;
    }
  }

  double sigcn = 1.0;

//---------------------------------------
//      Parameter Setting

  /*** Total excitation energy */
  prq.ex_total   = sys->ex_total;

  /*** Parameters for state density */
  for(int id=0 ; id<MAX_CHANNEL ; id++){
    if(!ncl[0].cdt[id].status) continue;
    int i = ncl[0].cdt[id].next;

    /*** Default single particle state density */
    spd[id].gz         = preqSingleStateDensity(ncl[i].za.getZ());
    spd[id].gn         = preqSingleStateDensity(ncl[i].za.getN());

    /*** Pairing energy */
    spd[id].pairing    = preqPairingDelta(ncl[i].za.getZ(),ncl[i].za.getA());

    /*** Potential well depth */
    spd[id].well_depth = preqPotentialDepth(sys->lab_energy,
                                            sys->incident.particle.getZ(),ncl[i].za.getA());
  }
  prq.spd = spd;


  /*** Initial number of particle/holes */
  int zp0, zh0, np0, nh0;
  if((sys->incident.particle.getZ() == 0) && (sys->incident.particle.getN() == 0)){
    /*** for photon-induced reaction, start with (1p,1h) in n-shell */
    zp0 = 0;  zh0 = 0;
    np0 = 1;  nh0 = 1;
  }
  else{
    zp0 = sys->incident.particle.getZ();   zh0 = 0;
    np0 = sys->incident.particle.getN();   nh0 = 0;
  }


//---------------------------------------
//      For Particle-Hole Configurations

  /*** Index i,j set to the number of holes */
  for(int i=0 ; i<MAX_PREEQ ; i++){
    exc.zp = zp0 + i;
    exc.zh = zh0 + i;
    for(int j=0 ; j<MAX_PREEQ ; j++){
      exc.np = np0 + j;
      exc.nh = nh0 + j;

      /*** Total exciton state density omega(p,h,Ex) */
      prq.omega_total = preqStateDensity(prq.ex_total,&prq.spd[0],&exc);

      if(prq.omega_total==0.0) continue;

      /*** Particle emission rate */
      wk[i][j] = preqEmissionRate(ncl,pdt,tc,&prq,&exc,wkc);

      /*** Effective squared matrix elements */
      int n = exc.zp + exc.zh + exc.np + exc.nh;
      preqMSquare(sys->target.getA()+sys->incident.particle.getA(),n,&prq);

      /*** Lambdas for all configurations */
      preqLifeTime(i,j,wk,tzp,tnp,tz0,tn0,tau,tap,&prq,&exc);
    }
  }


//---------------------------------------
//      Ocupation Probability for Each Configuration

  preqPOccupation(pop,tau,tap,tzp,tnp,tz0,tn0);

//---------------------------------------
//      Add to Population

  for(int i=0 ; i<MAX_PREEQ ; i++){
    exc.zp = zp0 + i;
    exc.zh = zh0 + i;
    for(int j=0 ; j<MAX_PREEQ ; j++){
      exc.np = np0 + j;
      exc.nh = nh0 + j;
      prq.omega_total = preqStateDensity(prq.ex_total,&prq.spd[0],&exc);
      if(prq.omega_total==0.0) continue;

      preqEmissionRate(ncl,pdt,tc,&prq,&exc,wkc);

      double taup = tau[i][j] * pop[i][j];

      for(int id=1 ; id<MAX_CHANNEL ; id++){
        if(!ncl[0].cdt[id].status) continue;
        int ip = ncl[0].cdt[id].next;
        for(int k=1 ; k<ncl[ip].ntotal ; k++){
          spc[id][k] += sigcn * wkc[id][k] * taup;
        }
      }
    }
  }


//---------------------------------------
//      Free Allocated Memory

  for(int j=0 ; j<MAX_CHANNEL ; j++){
    delete [] wkc[j];
  }
  for(int j=0 ; j<MAX_PREEQ   ; j++){
    delete [] pop[j];
    delete [] tzp[j];
    delete [] tnp[j];
    delete [] tz0[j];
    delete [] tn0[j];
    delete [] wk[j];
    delete [] tau[j];
    delete [] tap[j];
  }
}


/**********************************************************/
/*      Particle Emission Rate                            */
/*      calculated in the unit of h-bar [MeV sec]         */
/**********************************************************/
double preqEmissionRate(Nucleus *ncl, Pdata *pdt, Transmission **tc, Preeq *prq, Exconf *exc, double **wkc)
{
  double c  = AMUNIT/(HBARSQ*VLIGHTSQ*NORM_FACT*PI*PI);
  double wk = 0.0;

  for(int id=1 ; id<MAX_CHANNEL ; id++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) wkc[id][k] = 0.0;
    if(!ncl[0].cdt[id].status) continue;

    int j = ncl[0].cdt[id].next;

    Exconf res(exc->zp - pdt[id].particle.getZ(),exc->zh,
               exc->np - pdt[id].particle.getN(),exc->nh);
    if(res.zp<0 || res.np<0) continue;

    double rm = ncl[j].mass * pdt[id].mass
              /(ncl[j].mass + pdt[id].mass);

    double x  = c * rm * (2.0*pdt[id].spin+1.0)/prq->omega_total;

    for(int k=1 ; k<ncl[j].ntotal ; k++){
      double omega = preqStateDensity(ncl[j].excitation[k],&prq->spd[id],&res);
      wkc[id][k]   = x * tc[id][k].ecms * tc[id][k].sigr * omega;
      wk += wkc[id][k] * ncl[j].de;
    }
  }

  return(wk);
}


/**********************************************************/
/*      Lifetimes Tau of Exciton State                    */
/*      calculated in the unit of h-bar [MeV sec]         */
/**********************************************************/
void preqLifeTime(int i, int j, double **w,
                  double **z1, double **z2, double **z3, double **z4,
                  double **t1, double **t2, Preeq *prq, Exconf *exc)
{
  z1[i][j] = preqLambdaPlusZ(prq,exc);
  z2[i][j] = preqLambdaPlusN(prq,exc);
  z3[i][j] = preqLambdaZeroZ(prq,exc);
  z4[i][j] = preqLambdaZeroN(prq,exc);

  double d1 = z1[i][j] + z2[i][j]+ z3[i][j] + z4[i][j] + w[i][j];
  double d2 = z1[i][j] + z2[i][j]+ w[i][j];
  t1[i][j] = (d1>0.0) ? 1.0/d1 : 0.0;
  t2[i][j] = (d2>0.0) ? 1.0/d2 : 0.0;
}


/**********************************************************/
/*      Spin-Averaged Transmission Coefficients           */
/**********************************************************/
double preqTransmissionAverage(int spin2, int l, double *tj)
{
  double t=0.0;

  if(spin2==0)       t  = tj[0];
  else if(spin2==1)  t  = ( (l+1.0)*tj[0] + l *tj[1] )/(2.0*l+1.0);
  else if(spin2==2)  t  = ( (2.0*l+3.0)*tj[0]
                           +(2.0*l+1.0)*tj[1]
                           +(2.0*l-1.0)*tj[2] )/(6.0*l+3.0);
  return (t);
}

