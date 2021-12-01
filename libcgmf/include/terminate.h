/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file terminate.h

  \brief Abort/stop code execution

*/
#ifndef __TERMINATE_H__
#define __TERMINATE_H__

int     cgmTerminateCode      (std::string);
int     cgmTerminateCode      (std::string, int);
int     cgmTerminateCode      (std::string, double);

void    cgmDeleteAllocated    ();
void    cgmAllocateMemory     ();

#endif //__TERMINATE_H__
