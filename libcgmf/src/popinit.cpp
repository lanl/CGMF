/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file popinit.cpp

  \brief Store initial population into energy bin or level

*/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>

using namespace std;

#include "cgm.h"
#include "global.h"

/*** transmission coefficient index calculator */
static inline int tj_index(int a, int b, int c){
  return((a*3 + (1-b)*c)/2);
}

static void statPopinitLevden        (int, double, Nucleus *);
static void statPopinitSingle        (double, int, int, Nucleus *);
static void statPopinitTransmission  (int, Nucleus *, Nucleus *, Transmission *);
static void statWriteInitPopoulation (int, Parity *);

/**********************************************************/
/*      Store Initial Population at the Given Excitation  */
/**********************************************************/
void  statSetupInitialPopulationSpec(double j0, int p0, int k0, double sf, Transmission *tin)
{
  ncl[0].jmax=MAX_J-1;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++) ncl[0].pop[k][j].even = ncl[0].pop[k][j].odd  = 0.0;
  }

  if(ctl.init_pop==LEVDEN){
    statPopinitLevden(k0,sf,&ncl[0]);
  }else if(ctl.init_pop==TRANSMISSION){
    statPopinitTransmission(k0,&ncl[0],&ncl[1],tin);
  }else{
    statPopinitSingle(j0,p0,k0,&ncl[0]);
  }

  if(ctl.print_init_pop) statWriteInitPopoulation(ncl[0].jmax,ncl[0].pop[k0]);
}


/**********************************************************/
/*      Store Initial Population for Beta-Decay Case      */
/**********************************************************/
void  statSetupInitialPopulationBeta(double j0, int p0, int kmax)
{
  ncl[0].jmax=MAX_J-1;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++) ncl[0].pop[k][j].even = ncl[0].pop[k][j].odd  = 0.0;
  }

  double m = (j0 <= 0.5) ? 2.0 : 3.0;

  for(int k=0 ; k<kmax ; k++){
    if(j0 <= 0.5){
      statPopinitSingle(j0    ,p0,k,&ncl[0]);
      statPopinitSingle(j0+1.0,p0,k,&ncl[0]);
    }
    else{
      statPopinitSingle(j0-1.0,p0,k,&ncl[0]);
      statPopinitSingle(j0    ,p0,k,&ncl[0]);
      statPopinitSingle(j0+1.0,p0,k,&ncl[0]);
    }
    for(int j=0 ; j<=ncl[0].jmax ; j++){
      ncl[0].pop[k][j].even /= m;
      ncl[0].pop[k][j].odd  /= m;
    }
    if(ctl.print_init_pop) statWriteInitPopoulation(ncl[0].jmax, ncl[0].pop[k]);
  }
}


/**********************************************************/
/*      Output Initial Popolation                         */
/**********************************************************/
void statWriteInitPopoulation(int jmax, Parity *p)
{
  double sum = 0.0;
  for(int j=0 ; j<=jmax ; j++){
    cout << setw(11) << j;
    cout << setprecision(4) << setiosflags(ios::scientific)
         << setw(12) << p[j].even
         << setw(12) << p[j].odd << endl;
    sum += p[j].even + p[j].odd;
  }
  cout  << "# total    " << setw(12) << sum << endl;
}



/**********************************************************/
/*      Population Distributed Over the Spin States       */
/**********************************************************/
void statPopinitLevden(int k0, double sf, Nucleus *n)
{
  int jmax = 0;
  double sd, se=0.0, so=0.0;

  ldLevelDensity(n->excitation[0],(double)n->za.getA(),&n->ldp);

  /*** scale spin-cutoff parameter by input */
  if(sf > 0.0) n->ldp.spin_cutoff *= sf;
  //  cout << n->za.getA() << "  " << n->ldp.spin_cutoff << endl;

  for(int j=0 ; j<MAX_J ; j++){
    sd = ldSpinDistribution(j+halfint(n->lev[0].spin),n->ldp.spin_cutoff,n->ldp.sigma0);

    se += (n->pop[k0][j].even = sd * ldParityDistribution( 1));
    so += (n->pop[k0][j].odd  = sd * ldParityDistribution(-1));
    if(2*sd<SPIN_CUTOFF && j>2){
      jmax=j;
      break;
    }
  }
  sd = se+so;
  se = so= 0.0;
  for(int j=0;j<=jmax;j++){
    n->pop[k0][j].even /= sd;
    n->pop[k0][j].odd  /= sd;
    se += n->pop[k0][j].even;
    so += n->pop[k0][j].odd ;
  }
}


/**********************************************************/
/*      Popolation From Transmission Coefficients         */
/**********************************************************/
void statPopinitTransmission(int k0, Nucleus *n0, Nucleus *n1, Transmission *tin)
{
  int ig    = (int)(2.0*halfint(n0->lev[0].spin)); // half-int ofset for CN
  int it    = (int)(2.0*n1->lev[0].spin);          // target state spin
  int pt    = n1->lev[0].parity;                   // target state parity

  /*** Loop over CN J and Parity*/
  for(int j0=ig ; j0<=n0->jmax*2+ig ; j0+=2){
    int   jdx = (j0-ig)/2;
    for(int p0=-1 ; p0<=1 ; p0+=2){

      /*** For entrance channel Tlj coupled to JP*/
      for(int lp0=0 ; lp0<=tin->lmax*2 ; lp0+=2){
        if(parity(lp0) != p0*pt) continue;

        for(int sp0=1 ; sp0>=-1 ; sp0-=2){
          int jp0 = lp0 + sp0;
          if(jp0<0) continue;

          if( abs(jp0-it)>j0 || j0>(jp0+it) ) continue;

          double tlj  = tin->tran[tj_index(lp0,sp0,1)];
          if(tlj==0.0) continue;

          if(p0==1)  n0->pop[k0][jdx].even += (j0+1.0)*tlj;
          else       n0->pop[k0][jdx].odd  += (j0+1.0)*tlj;
        }
      }
    }
  }

  double sd = 0.0;
  for(int j=0 ; j<=n0->jmax ; j++){
    sd += n0->pop[k0][j].even;
    sd += n0->pop[k0][j].odd ;
  }
  for(int j=0 ; j<=n0->jmax ; j++){
    n0->pop[k0][j].even /= sd;
    n0->pop[k0][j].odd  /= sd;
  }
}


/**********************************************************/
/*      Popolation Localized on a Single JP State         */
/**********************************************************/
void statPopinitSingle(double j0, int p0, int k0, Nucleus *n)
{
  for(int j=0 ; j<MAX_J ; j++){
    if(j == (int)j0){
      if(p0 < 0){
        n->pop[k0][j].even  = 0.0;
        n->pop[k0][j].odd  += 1.0;
      }else{
        n->pop[k0][j].even += 1.0;
        n->pop[k0][j].odd   = 0.0;
      }
    }
  }
}


/***********************************************************/
/*     Find Nearest Discrete Level, and  Store Population  */
/***********************************************************/
int  statStoreLevelExcite(Nucleus *n)
{
  double d=100.0,x;
  int k=0;

  if( (n->ndisc == 1) && (n->ncont == 0) ) return(k);

  for(int i=0 ; i<n->ndisc-1 ; i++){
    x = fabs(n->max_energy - n->lev[i].energy);
    if(d>x){
      k = i;
      d = x;
    }
  }
  if(k>0) n->lpop[k] = 1.0;

  if(ctl.print_init_pop){
      cout << setw(11) << k;
      cout << setprecision(4) << setiosflags(ios::scientific)
           << setw(12) << n->lev[k].energy
           << setw(12) << n->lpop[k] << endl;
  }

  return(k+1);
}


/**********************************************************/
/*      Fill Zeros in Population Array                    */
/**********************************************************/
void statClearPopulation(Nucleus *n)
{
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++) {
      n->pop[k][j].even = n->pop[k][j].odd = 0.0;
    }
  }
  for(int k=0 ; k<MAX_LEVELS ; k++) n->lpop[k] = 0.0;
}
