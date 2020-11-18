/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file evapinterface.h

  \brief Interface for evaporation spectrum

*/

#ifndef __EVAPINTERFACE_H__
#define __EVAPINTERFACE_H__

#include "structur.h"

void evapInterface(Pdata *pdt,Transmission *tc,int c0,int k0,double *sp);

#endif //__EVAPINTERFACE_H__
