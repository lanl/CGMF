/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file masstable.cpp

  \brief Nuclear masses

*/

/******************************************************************************/
/*        Find nuclear mass from Mass Table. Allows for easy interchange      */
/*        of the theoretical mass table. Ensures Sn or Q-values are only      */
/*        determined from a difference between theory-theory or exp-exp.      */
/******************************************************************************/

/* C INCLUDES */
#include <iostream>
#include <sstream>

/* HEADER INCLUDES */
#include "masstable.h"
#include "terminate.h"

//#include "masstable_aw95.h"      // Audi Wapstra 1995 table
//#include "masstable_ripl2.h"     // AW95 + FRDM95 from RIPL2
//#include "masstable_ripl3.h"     // AW03 + FRDM95 from RIPL3
#include "masstable_audi2011.h"    // AW11 + FRDM95 from RIPL3

using namespace std;

/*******************************************************************************
 * mass_excess
 *------------------------------------------------------------------------------
 * Calculates the mass excess for a given Z and A.
 ******************************************************************************/
double mass_excess(int z, int a) {

    double mx = 0.0;
    unsigned int za = z*1000 + a;

    bool found = false;
    // Search for this Z and A in the mass table
    for(int i=0;i<nMassTable;i++){
	if(MassTable[i].za == za){
	    found = true;
	    mx = MassTable[i].mass; // Pull the mass-excess
	    break;
	}
    }

    // Send message if Z,A wasn't found
    if(!found){
	ostringstream os;
	os << "mass data for Z " << z << " - A " << a << " not found";
	cgmTerminateCode(os.str());
    }

    return(mx);
}

