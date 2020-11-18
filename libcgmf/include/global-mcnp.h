/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file global-mcnp.h

  \brief Global variables to share with MCNP

*/

#ifndef __GLOBAL_MCNP_H__
#define __GLOBAL_MCNP_H__

#include "structur.h"

extern double **pfnInfo; // << 1.0.6 >>

//extern	FissionEvents		*fissionEvents;
extern	Nucleus					*ncl;                 // 0: parent, 1 - MAX_COMPOUND-1: daughters etc


#endif //__GLOBAL_MCNP_H__
