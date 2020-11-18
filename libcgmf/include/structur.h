/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file structur.h

  \brief Class to define the structure of a nucleus

*/

#ifndef __STRUCTUR_H__
#define __STRUCTUR_H__

/*
   structure.h :
        common structure definition
*/

#include "dimensions.h"
#include "gdr.h"
/*************************************************/
/*                     Enum                      */
/*************************************************/

enum particle{
  gammaray =0, neutron =1, proton  =2, alpha   =3, deuteron=4,
  triton   =5, helion  =6, fission =7, unknown =8};
typedef enum particle Particle;


typedef enum {sumall=0, hauser=1} Statcalcmode;


/*************************************************/
/*                  Class                        */
/*************************************************/

/*!

 \brief A very complicated (and slow?) way to define a nucleus

 \todo Add complete documentation.
 
 */
class ZAnumber{
 private:
    unsigned int Z;
    unsigned int A;
 public:
    ZAnumber(){
      Z = 0;
      A = 0;
    }
    ZAnumber(int z, int a){
      Z = z;
      A = a;
    }
    void setZA(int z, int a){
      Z = z;
      A = a;
    }
    unsigned int getZ(){ return (Z); }
    unsigned int getA(){ return (A); }
    unsigned int getN(){ return (A-Z); }
    ZAnumber operator+(ZAnumber x){
      ZAnumber y;
      y.Z = Z + x.Z;
      y.A = A + x.A;
      return y;
    }
    ZAnumber operator-(ZAnumber x){
      ZAnumber y;
      y.Z = Z - x.Z;
      y.A = A - x.A;
      return y;
    }
};


/****************************/
/*   Array for Parites      */
/****************************/
class Parity{
 public:
    double even              ;
    double odd               ;

    Parity(){
      even = 0.0;
      odd  = 0.0;
    }
};

/****************************/
/*   Particle Data          */
/****************************/
class Pdata{
 public:
    Particle  particleID     ;     /* particle identifier              */
    ZAnumber  particle       ;     /* particle mass and atomic number  */
    unsigned int omp         ;     /* optical potential index          */
    double    spin           ;     /* spin                             */
    double    mass           ;     /* exact particle mass              */
    double    mass_excess    ;     /* particle mass excess             */

    Pdata(){
      omp         = 0;
      spin        = 0.0;
      mass        = 0.0;
      mass_excess = 0.0;
    }
};

/****************************/
/*    Discrete Levels       */
/****************************/
class Level{
 public:
    double     energy        ;     /* Level energy                     */
    double     spin          ;     /* Level spin                       */
    int        parity        ;     /* Parity                           */
    int        ngamma        ;     /* Number of Gamma-rays             */
    int       *fstate        ;     /* Final state index for g-decay    */
    double    *branch        ;     /* Branching ratio                  */
    double     halflife      ;     /* half-life of the state           */
    double    *gratio        ;     /* net gamma-ray considering ICC    */

    Level(){
      energy   = 0.0;
      spin     = 0.0;
      parity   = 0;
      ngamma   = 0;
      halflife = 0.0;
      fstate   = new int    [MAX_GAMMA_BRANCH];
      branch   = new double [MAX_GAMMA_BRANCH];
      gratio   = new double [MAX_GAMMA_BRANCH];
    }
    Level(const Level& that){
      gratio   = new double [MAX_GAMMA_BRANCH];
      branch   = new double [MAX_GAMMA_BRANCH];
      fstate   = new int    [MAX_GAMMA_BRANCH];
      for(int i=0 ; i<MAX_GAMMA_BRANCH ; i++){
        gratio[i] = that.gratio[i];
        branch[i] = that.branch[i];
        fstate[i] = that.fstate[i];
      }
      halflife = that.halflife;
      ngamma   = that.ngamma;
      parity   = that.parity;
      spin     = that.spin;
      energy   = that.energy;
    }
    Level& operator=(const Level& that){
      double* local_gratio = new double [MAX_GAMMA_BRANCH];
      double* local_branch = new double [MAX_GAMMA_BRANCH];
      int*    local_fstate = new int    [MAX_GAMMA_BRANCH];
      for(int i=0 ; i<MAX_GAMMA_BRANCH ; i++){
        local_gratio[i] = that.gratio[i];
        local_branch[i] = that.branch[i];
        local_fstate[i] = that.fstate[i];
      }
      delete [] gratio;
      delete [] branch;
      delete [] fstate;
      gratio = local_gratio;
      branch = local_branch;
      fstate = local_fstate;
      halflife = that.halflife;
      ngamma   = that.ngamma;
      parity   = that.parity;
      spin     = that.spin;
      energy   = that.energy;
      return *this;
    }
    ~Level(){
      delete [] gratio;
      delete [] branch;
      delete [] fstate;
    }
};

/****************************/
/*   Channel Depend  Data   */
/****************************/
class Channel{
 public:
    Particle particleID      ;     /* particle identifier              */
    int      spin2           ;     /* particle spin x2                 */
    int      next            ;     /* index of residual nucleus        */
    int      lmax            ;     /* maximal angular momentum         */
    double   energy          ;     /* energy of the emitted particle   */
    double   excitation      ;     /* energy of the residual nucleus   */
    double   binding_energy  ;     /* binding energy of the particle   */
    bool     status          ;     /* channel open / close             */

    Channel(){
      particleID     = unknown;
      spin2          = 0;
      next           = 0;
      lmax           = 0;
      energy         = 0.0;
      excitation     = 0.0;
      binding_energy = 0.0;
      status         = false;
    }
};


/****************************/
/*    Level Density         */
/****************************/
class LevelDensity{
 public:
    double     match_energy  ;     /* Matching energy                  */
    double     shell_correct ;     /* Shell correction energy          */
    double     temperature   ;     /* T parameter                      */
    double     E0            ;     /* Eergy shift                      */
    double     a             ;     /* a parameter                      */
    double     pairing_energy;     /* Pairing energy                   */
    double     spin_cutoff   ;     /* Spin Cut-off factor              */
    double     sigma0        ;     /* spin Cut-off factor from levels  */

    LevelDensity(){
      match_energy   = 0.0;
      shell_correct  = 0.0;
      temperature    = 0.0;
      E0             = 0.0;
      a              = 0.0;
      pairing_energy = 0.0;
      spin_cutoff    = 0.0;
      sigma0         = 0.0;
    }
};

/****************************/
/*    Nucleus Data          */
/****************************/
class Nucleus{
 public:
    ZAnumber   za            ;     /* Z and A numbers                  */
    double     max_energy    ;     /* initial excitation energy        */
    double     mass          ;     /* exact mass                       */
    double     mass_excess   ;     /* mass excess                      */
    double     de            ;     /* energy bin width in continuum    */
    double     gstrength     ;     /* gamma-ray strength function      */
    LevelDensity ldp         ;     /* level density parameter          */
    Level     *lev           ;     /* discrete level data              */
    double    *binwidth      ;     /* energy bin width in continuum    */
    double    *excitation    ;     /* excitation energy                */
    Parity   **density       ;     /* level density                    */
    Parity   **pop           ;     /* population                       */
//    Parity   **dpop          ;     /* population increment             */
    double    *lpop          ;     /* leve population                  */
    Channel   *cdt           ;     /* channel data                     */
    int        jmax          ;     /* cut-off of angular momentum      */
    int        ndisc         ;     /* number discrete levels           */
    int        ncont         ;     /* number of continuum enerby bins  */
    int        ntotal        ;     /* total number of enerby bins      */

    Nucleus(){
      max_energy     = 0.0;
      mass           = 0.0;
      mass_excess    = 0.0;
      de             = 0.0;
      gstrength      = 0.0;
      jmax           = 0  ;
      ndisc          = 0  ;
      ncont          = 0  ;
      ntotal         = 0  ;

      density        = new Parity *[MAX_ENERGY_BIN];
      pop            = new Parity *[MAX_ENERGY_BIN];
      for(int j=0 ; j<MAX_ENERGY_BIN ; j++){
        density[j]     = new Parity[MAX_J];
        pop[j]         = new Parity[MAX_J];
      }
      lpop           = new double[MAX_LEVELS];
      binwidth       = new double[MAX_ENERGY_BIN];
      excitation     = new double[MAX_ENERGY_BIN];
      lev = new Level [MAX_LEVELS];
      cdt = new Channel [MAX_CHANNEL];
    }
    ~Nucleus(){
      for(int j=0 ; j<MAX_ENERGY_BIN ; j++){
        delete [] density[j];
        delete [] pop[j];
      }
      delete [] density;
      delete [] pop;

      delete [] lpop;
      delete [] binwidth;
      delete [] excitation;
      delete [] lev;
      delete [] cdt;
    }
};

/****************************/
/*   Cross Sections         */
/****************************/
class CrossSection{
 public:
    double  total            ;     /* total cross section              */
    double  elastic          ;     /* elastic cross section            */
    double  reaction         ;     /* non-elastic cross section        */

    CrossSection(){
      total     = 0.0;
      elastic   = 0.0;
      reaction  = 0.0;
    }
};


/**************************************/
/*      Transmission coefficients     */
/**************************************/
class Transmission{
 public:
    int       lmax           ;     /* max L                  */
    double    ecms           ;     /* cms energy             */
    double    sigr           ;     /* reaction cross section */
    double   *tran           ;     /* Tj [l+s][l][l-s]       */

    Transmission(){
      lmax = 0;
      ecms = 0.0;
      sigr = 0.0;
      tran = NULL;
    }
};


extern  Nucleus      *ncl;
extern  double       *spc[SPECTRA_OUTPUT];


#endif //__STRUCTUR_H__
