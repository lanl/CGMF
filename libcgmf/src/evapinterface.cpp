/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file evapinterface.cpp

  \brief Interface to compute an evaporation excitation spectrum

*/

#include <iostream>

using namespace std;

#include "cgm.h"
#include "evapinterface.h"
#include "terminate.h"

void evapInterface(Pdata *pdt,Transmission *tc,int c0,int k0,double *sp){
   
   /*** calculate inverse reaction cross sections */
//   cout << ncl[0].max_energy << endl;
   statStoreContinuumTransmission(c0,k0,pdt,tc);
//   statStoreDiscreteTransmission(c0,k0,pdt,td);
  
   fill_n (sp, MAX_ENERGY_BIN, 0.0);
   
   /*** Construct the neutron spectrum ****/
   double s=0.;
   for (int k=0; k<k0; ++k) {
      tc[k].ecms=0.;
   }
   double ecms0=0.;
   for (int k=k0; k<ncl[c0+1].ntotal; ++k) {
      double ecms=tc[k].ecms;
      if(ecms>0){
         sp[k]=tc[k].sigr*ecms*ldLevelDensity(ncl[c0+1].excitation[k], (double)ncl[c0+1].za.getA(), &ncl[c0+1].ldp);
         s+=sp[k];
      }
      ecms0=ecms;
   }
   
   for (int k=k0; k<ncl[c0+1].ntotal; ++k) { // normalize the spectrum to one
     sp[k]/=s;
     //     if(tc[k].ecms>0.)
     //       cout <<  tc[k].ecms << " " << sp[k] << endl;
   }

   return;
}
