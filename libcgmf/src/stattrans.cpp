/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file stattrans.cpp

  \brief Set up transmission coefficients

*/

#include <iostream>
#include <cmath>

using namespace std;

#include "physics.h"
#include "cgm.h"
#include "masstable.h"


/**********************************************************/
/*      Store Particle Continuum Transmission into Array  */
/**********************************************************/
void statStoreContinuumTransmission(int ip, int kp, Pdata *p, Transmission *tc)
{
  CrossSection cx;
  double mu = 0.0;
  /*** Clear transmission array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    tc[k].lmax = 0;
    tc[k].ecms = 0.0;
//    for(int l=0 ; l<3*MAX_J ; l++) tc[k].tran[l] = 0.0;
  }

  double sp = ncl[ip].cdt[neutron].binding_energy;
  int    id = ncl[ip].cdt[neutron].next;
	
  if( id != -1 ) {
    mu = ncl[id].mass * p->mass / (ncl[id].mass + p->mass);
	
    for(int k=0 ; k<ncl[id].ntotal ; k++){
      tc[k].ecms = ncl[ip].excitation[kp] - ncl[id].excitation[k] - sp;
      if(tc[k].ecms<0.0) continue;

      tc[k].lmax = omCalc(tc[k].ecms, p, &ncl[id].za, mu, tc[k].tran, &cx);
      tc[k].sigr = cx.reaction;
    }
  }
}


/**********************************************************/
/*      Store Particle Discrete Transmission into Array   */
/**********************************************************/
void statStoreDiscreteTransmission(int ip, int kp, Pdata *p, Transmission *td)
{
  CrossSection cx;
  double mu = 0.0;
  /*** Clear transmission array */
  for(int k=0 ; k<MAX_LEVELS ; k++){
    td[k].lmax = 0;
    td[k].ecms = 0.0;
    // for(int l=0 ; l<3*MAX_J ; l++) td[k].tran[l] = 0.0;
  }

  double sp = ncl[ip].cdt[neutron].binding_energy;
  int    id = ncl[ip].cdt[neutron].next;
  if( id != -1 ) {
    mu = ncl[id].mass * p->mass / (ncl[id].mass + p->mass);

    /*** Transition to discrete levels */
    for(int k=0 ; k<ncl[id].ndisc ; k++){
      td[k].ecms = ncl[ip].excitation[kp] - ncl[id].lev[k].energy - sp;
      if(td[k].ecms<0.0) continue;
      td[k].lmax = omCalc(td[k].ecms, p, &ncl[id].za, mu, td[k].tran, &cx);
      td[k].sigr = cx.reaction ;
    }
  }
}


/**********************************************************/
/*  Store Gamma-ray Transmissions into 2-dim Array        */
/**********************************************************/
void statStoreGammaTransmission(int k0, double **tg, Nucleus *n)
{
  /*** Clear transmission array
  for(int i=0 ; i<MAX_MULTIPOL ; i++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tg[i][k] = 0.0;
  }*/

  /*** Neutron separation energy for temperature dependent Gamma */
  double sn = n->cdt[neutron].binding_energy;
  if(sn==0.0){
    double mx = mass_excess(n->za.getZ(),n->za.getA()-1);
    sn = mx - n->mass_excess + ENEUTRON;
  }

  for(int k1=k0+1 ; k1<n->ncont ; k1++){
    double eg = n->excitation[k0] - n->excitation[k1];
    double ex = ( eg<=sn ) ? sn-eg : 0.0;

    tg[E1][k1] = gdrGammaTransmission(GL,E1_MULTIPOL,eg,E1,n->ldp.a,ex);
    tg[M1][k1] = gdrGammaTransmission(SL,M1_MULTIPOL,eg,M1,0.0,0.0);
    tg[E2][k1] = gdrGammaTransmission(SL,E2_MULTIPOL,eg,E2,0.0,0.0);

/*
    double s = (tg[E1][k1] + tg[M1][k1])/(eg*eg*eg*PI2);
    cout <<k0 <<"  " << k1 << "  " << eg << "  " << ex << "  " << s << endl;
*/
  }
}
