/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file masstable.h

  \brief Mass excess table

*/

#ifndef __MASSTABLE_H__
#define __MASSTABLE_H__

/* C INCLUDES */
#include <string>
#include <vector>

// Default mass excess table (hardcoded in masstable_audi2011.h)
class MassExcess{
    public:
	unsigned int za;    // Z*1000 + A
	float        mass;  // mass excess
};

double  mass_excess(int,int);

#endif //__MASSTABLE_H__
