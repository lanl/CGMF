/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file specevap.cpp

  \brief Particle emission using evaporation model

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "terminate.h"

static void evapConstantTemperature(int, double, double, Transmission *, Nucleus *, double *);
static void evapLevelDensity(int, double, double, Transmission *, Nucleus *, double *);

const string pmodel = "levelden";
//const string pmodel = "constemp";
//const string pmodel = "fermigas";


/**********************************************************/
/*      Neutron Emission Spectrum Using Evaporation       */
/**********************************************************/
void specEvaporation(int nc, Pdata *pdt, double **spc)
{
  Transmission  *tc = NULL;
  double        *gn[2];
  
  //---------------------------------------
  //      Allocate Transmission Data
  
  try{
    tc = new Transmission [MAX_ENERGY_BIN];
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
    for(int i=0 ; i<2 ; i++)              gn[i]      = new double [MAX_ENERGY_BIN];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in spectraEvaporation");
  }
  
  /*** clear total spectrum array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[0][k] = spc[1][k] = spc[2][k] = spc[3][k] = 0.0;
  
  /*** initial population stored in J=0 */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++) ncl[0].pop[k][j].even = ncl[0].pop[k][j].odd  = 0.0;
  }
  ncl[0].pop[0][0].even = 1.0;
  ncl[0].pop[0][0].odd  = 0.0;
  
  
  //---------------------------------------
  //      Nucleus Decay Probability by Weisskopf
  
  /*** for all CN */
  for(int c0=0; c0<nc ; c0++){
    
    /*** calculate transmission coefficients for neutron channel */
    if(!ncl[c0].cdt[neutron].status) continue;
    
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) gn[0][k] = gn[1][k] = 0.0;
    
    statStoreContinuumTransmission(c0,0,&pdt[neutron],tc);
    
    int c1 = ncl[c0].cdt[neutron].next;
    
    /*** loop over parent excitation */
    for(int k0=0 ; k0<ncl[c0].ntotal ; k0++){
      
      double ex = ncl[c0].excitation[k0] - ncl[c0].cdt[neutron].binding_energy;
      if(ex <= 0.0) continue;
      
      double p0 = ncl[c0].pop[k0][0].even;
      if(p0 == 0.0) continue;
      
      if(pmodel == "levelden") evapLevelDensity(k0,p0,ex,tc,&ncl[c1],gn[1]);
      else                     evapConstantTemperature(k0,p0,ex,tc,&ncl[c1],gn[1]);
      
    }
    /*** cumulative neutron spectra */
    for(int k=0 ; k<ncl[c0].ntotal ; k++) spc[1][k] += gn[1][k];
    
    if(ctl.print_indivspec) cgmPrintSpectra(false,ncl[0].de,gn);
  }
  
  //---------------------------------------
  //      Free Allocated Memory
  
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[k].tran;
  for(int i=0 ; i<2 ; i++)              delete [] gn[i];
  delete [] tc;
}



/**********************************************************/
/*      Simple Evaporation From Fixed State               */
/**********************************************************/
void specEvaporationBin(int c0, Pdata *pdt, double **spc)
{
  if(!ncl[c0].cdt[neutron].status) return;
  
  Transmission  *tc = NULL;
  
  //---------------------------------------
  //      Allocate Transmission Data
  
  try{
    tc = new Transmission [MAX_ENERGY_BIN];
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in spectraEvaporation");
  }
  
  /*** clear total spectrum array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[0][k] = spc[1][k] = 0.0;
  
  
  //---------------------------------------
  //      Nucleus Decay Probability by Weisskopf
  
  
  /*** calculate transmission coefficients for neutron channel */
  statStoreContinuumTransmission(c0,0,&pdt[neutron],tc);
  
  int c1 = ncl[c0].cdt[neutron].next;
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int j=0 ; j<MAX_J ; j++){
      ncl[c0].pop[k][j].even = ncl[c0].pop[k][j].odd  = 0.0;
      ncl[c1].pop[k][j].even = ncl[c1].pop[k][j].odd  = 0.0;
    }
  }
  
  double ex = ncl[c0].excitation[0] - ncl[c0].cdt[neutron].binding_energy;
  if(ex > 0.0){
    if(pmodel == "levelden") evapLevelDensity(0,1.0,ex,tc,&ncl[c1],spc[1]);
    else                     evapConstantTemperature(0,1.0,ex,tc,&ncl[c1],spc[1]);
  }
  
  //---------------------------------------
  //      Free Allocated Memory
  
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[k].tran;
  delete [] tc;
}


/**********************************************************/
/*      Standard Evaporation Spectrum with Constant T     */
/**********************************************************/
void evapConstantTemperature(int k0, double p0, double ex,
                             Transmission *tc, Nucleus *n1, double *sn)
{
  double t  = n1->ldp.temperature;
  
  if(pmodel == "fermigas"){
    double a  = ldDensityParameter(ex, (double)n1->za.getA(), &n1->ldp);
    t = sqrt(ex / a);
  }
  
  double st = 0.0;
  for(int k1=k0 ; k1<n1->ntotal ; k1++){
    int    kp = k1-k0;
    if(tc[kp].lmax <= 0) continue;
    
    double ep = ex - n1->excitation[k1];
    st += ep*exp(-ep/t) * tc[kp].sigr;
  }
  if(st == 0.0) return;
  
  for(int k1=k0 ; k1<n1->ntotal ; k1++){
    int    kp = k1-k0;
    if(tc[kp].lmax <= 0) continue;
    
    double ep = ex - n1->excitation[k1];
    double ds = p0 * ep * exp(-ep/t) * tc[kp].sigr / st;
    
    sn[kp] += ds;
    n1->pop[k1][0].even += ds;
  }  
}


/**********************************************************/
/*      Evaporation Spectrum based on Level Density       */
/**********************************************************/
void evapLevelDensity(int k0, double p0, double ex,
                      Transmission *tc, Nucleus *n1, double *sn)
{
  double st = 0.0;
  for(int k1=k0 ; k1<n1->ntotal ; k1++){
    int    kp = k1-k0;
    if(tc[kp].lmax <= 0) continue;
    
    double ep = ex - n1->excitation[k1];
    double r  = ldLevelDensity(n1->excitation[k1],(double)n1->za.getA(),&n1->ldp);
    st += ep* r * tc[kp].sigr;
  }
  if(st == 0.0) return;
  
  for(int k1=k0 ; k1<n1->ntotal ; k1++){
    int    kp = k1-k0;
    if(tc[kp].lmax <= 0) continue;
    
    double ep = ex - n1->excitation[k1];
    double r  = ldLevelDensity(n1->excitation[k1],(double)n1->za.getA(),&n1->ldp);
    double ds = p0 * ep * r * tc[kp].sigr / st;
    
    sn[kp] += ds;
    n1->pop[k1][0].even += ds;
  }  
}

