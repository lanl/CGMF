/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file gamcas.cpp

  \brief Gamma-ray cascading

*/

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "config.h"


/*************************************************/
/*  Gamma Cascading From Each Level              */
/*************************************************/
int specGammaCascade(int nstart, double *spc, Nucleus *n)
{
  double s=0.0;
  int    k,m=0;

  for(int i0=nstart-1 ; i0>0 ; i0--){
    for(int j=0 ; j<n->lev[i0].ngamma ; j++){
      int i1 = n->lev[i0].fstate[j];
      n->lpop[i1] += (s=n->lpop[i0]*n->lev[i0].branch[j]);
      if(s>0.0){
        double e = n->lev[i0].energy - n->lev[i1].energy;
        if( (k = specFindEnergyBin(e,n->de))>=0 ){
          if(INCLUDE_INTERNAL_CONVERSION) s *= n->lev[i0].gratio[j];
          spc[k] += s;
          m++;
        }
      }
    }
  }
  return (m);
}


/*************************************************/
/*  Gamma Branching From a Given Level           */
/*************************************************/
int specGammaDiscrete(int nstart, double *spc, Nucleus *n)
{
  double s=0.0;
  int    k,m=0;
  int    i0=nstart-1;

  for(int j=0 ; j<n->lev[i0].ngamma ; j++){
    int i1 = n->lev[i0].fstate[j];
    n->lpop[i1] += (s=n->lpop[i0]*n->lev[i0].branch[j]);
    if(s>0.0){
      double e = n->lev[i0].energy - n->lev[i1].energy;
#ifdef HAVE_PRIVATE_ENERGY_GRID
      k = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,e,n->binwidth);
#else
      k = specFindEnergyBin(e,n->de);
#endif
      if(k>=0){
        if(INCLUDE_INTERNAL_CONVERSION) s *= n->lev[i0].gratio[j];
        spc[k] += s;
        m++;
      }
    }
  }
  return (m);
}
