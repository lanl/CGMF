/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file global.h

  \brief Global variables and routines

*/
#ifndef __GLOBAL_H__
#define __GLOBAL_H__

#define OUT_SPEC       0x0001    //  1
#define OUT_INIPOP     0x0002    //  2
#define OUT_HIST       0x0004    //  4
#define OUT_BROADENED  0x0008    //  8
#define OUT_INDIVSPEC  0x0010    // 16
//#define              0x0020    // 32
//#define             0x0040
//#define             0x0080
#define CALC_MC     0x0100
#define CALC_TRAN   0x0400

typedef enum {SINGLE, LEVDEN, BETA, TRANSMISSION} InitialPop;

#define DEFAULT_CAL (OUT_SPEC)


class Calculation{ 
 public:
    bool print_history    ;
    bool print_spectrum   ;
    bool print_init_pop   ;
    bool print_broadened  ;
    bool print_indivspec  ;
    bool calc_montecarlo  ;
    bool calc_entrance    ;
    InitialPop init_pop   ;
};

extern Calculation ctl;

extern int mcl_nlines;
extern double mcl_glines[];

#endif //__GLOBAL_H__
