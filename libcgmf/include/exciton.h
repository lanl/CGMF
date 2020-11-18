/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file exciton.h

  \brief Prototype of functions for exciton model calculations

*/


#ifndef __EXCITON_H__
#define __EXCITON_H__


#define MAX_PREEQ  6

/**************************************/
/*      Exciton Configuration         */
/**************************************/
class Exconf{
 public:
    int        zp            ;     /* proton particle number              */
    int        zh            ;     /* proton hole number                  */
    int        np            ;     /* neutron particle number             */
    int        nh            ;     /* neutron hole number                 */

    Exconf(){
      zp            =0   ;
      zh            =0   ;
      np            =0   ;
      nh            =0   ;
    }
    Exconf(int p1, int h1, int p2, int h2){
      zp            =p1  ;
      zh            =h1  ;
      np            =p2  ;
      nh            =h2  ;
    }
    void setExconf(int p1, int h1, int p2, int h2){
      zp            =p1  ;
      zh            =h1  ;
      np            =p2  ;
      nh            =h2  ;
    }
};


/**************************************/
/*      Singple Particle Density      */
/**************************************/
class SPdens{
 public:
    double     gz            ;     /* state density parameter for proton  */
    double     gn            ;     /* state density parameter for neutron */
    double     pairing       ;     /* pairing energy                      */
    double     well_depth    ;     /* V, finite potential correction      */

    SPdens(){
      gz            =0.0 ;
      gn            =0.0 ;
      pairing       =0.0 ;
      well_depth    =0.0 ;
    }
};


/**************************************/
/*      Preequilibrium Systyem Data   */
/**************************************/
class Preeq{
 public:
    double     ex_total      ;     /* maximum excitation energy           */
    double     omega_total   ;     /* total state density for p-h         */
    double     m2zz          ;     /* averaged matrix elements for p-p    */
    double     m2nn          ;     /* averaged matrix elements for n-n    */
    double     m2zn          ;     /* averaged matrix elements for p-n    */
    double     m2nz          ;     /* averaged matrix elements for n-p    */
    SPdens     *spd          ;     /* pointer to s-p density array        */

    Preeq(){
      ex_total      =0.0 ;
      omega_total   =0.0 ;
      m2zz          =0.0 ;
      m2nn          =0.0 ;
      m2zn          =0.0 ;
      m2nz          =0.0 ;
    }
};


/**************************************/
/*      exciparm.cpp                  */
/**************************************/

double  preqSingleStateDensity     (int);
double  preqPairingDelta           (int, int);
double  preqPairingCorrection      (int, double, SPdens *);
double  preqPauliCorrection        (int, int, int, int, double, double);
double  preqMSquare                (int, int, Preeq *);
double  preqFiniteWell             (int, int, double, double);
double  preqPotentialDepth         (double, int, int);
double  preqStateDensity           (double, SPdens *, Exconf *);


/**************************************/
/*      excilambda.cpp                */
/**************************************/

double  preqLambdaPlusZ            (Preeq *, Exconf *);
double  preqLambdaPlusN            (Preeq *, Exconf *);
double  preqLambdaZeroZ            (Preeq *, Exconf *);
double  preqLambdaZeroN            (Preeq *, Exconf *);


/**************************************/
/*      excipop.cpp                   */
/**************************************/

void    preqPOccupation            (double **,double **,double **,
                                    double **,double **,double **,double **);

#endif //__EXCITON_H__
