/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file physics.h

  \brief Physics constants

*/

#ifndef __PHYSICS_H__
#define __PHYSICS_H__

// General physics constants

const double pi                           = 3.14159265359;
const double twopi                        = 6.28318530718;   // 2*pi
const double elementaryCharge             = 1.602176462e-19; // [C]
const double inverseFineStructureConstant = 137.03599976;
const double amu                          = 931.494013e+6;   // atomic mass unit [eV]
const double amuMeV                       = 931.494013;      // [MeV]
const double PlanckConstant               = 4.13566727e-15;  // [eV.s]
const double PlanckConstantOverTwoPi      = 6.58211889e-16;  // Planck constant / (2 pi) [eV.s]
const double BoltzmannConstant            = 8.617342e-5;     // [eV/K]
const double speedOfLight                 = 299792458.0;     // [m/s]
const double AvogadroNumber               = 6.02214199e23;   // [/mol]

const double hbarc = PlanckConstantOverTwoPi*speedOfLight*1e9;  // Hbar.c [MeV.fm]


/****************************/
/*      GENERAL             */
/****************************/

const double EXP1       =  2.71828182845904523536  ;  /* Napier's constant    */
const double PI         =  3.14159265358979323846  ;  /* circular constant    */
const double PI2        =  6.28318530717958647692  ;  /* PI*2                 */
const double PI4        = 12.56637061435917295384  ;  /* PI*4                 */
const double PI_2       =  1.57079632679489661923  ;  /* PI/2                 */
const double PI_4       =  0.785398163397448309616 ;  /* PI/4                 */
const double PI4_3      =  4.18879020478639053     ;  /* PI*4/3               */
const double SQRTPI     =  1.77245385090552        ;  /* sqrt(PI)             */
const double SQRT_1_PI  =  0.564189583547756286948 ;  /* sqrt(1/PI)           */
const double SQRT_2_PI  =  0.797884560802865406    ;  /* sqrt(2/PI)           */
const double LOG_2PI    =  1.83787706640934548     ;  /* log(2PI)             */
const double LOG_2      =  0.6931471805599453      ;  /* log(2.0)             */

//--- Masses [amu] ---

const double NEUTRONMASS  = 1.00866491578;
const double ELECTRONMASS = 5.48579911e-4;
const double PROTONMASS   = 1.00727646688;	
const double DEUTERONMASS = 2.01355321271;
const double TRITONMASS   = 3.016049268;
const double HE3MASS      = 3.01493223469;
const double ALPHAMASS    = 4.0015061747;

//--- Energies needed to break particles into their constituent nucleons [MeV] ---
const double DEUTERONBREAKUPENERGY = 2.22;
const double TRITONBREAKUPENERGY   = 8.48;
const double HE3BREAKUPENERGY      = 7.72;
const double ALPHABREAKUPENERGY    = 28.3;

//--- Conversion Factors ---
const double AMU2EV = 9.31494013e+8; //-- amu -> eV

/****************************/
/*      DEFAULT PARAMETER   */
/****************************/

// const double NORM_FACT      = 10.0 ; /* conversion factor to [mb]             */
// const double SPIN_CUTOFF    = 1e-8 ; /* cut-off of total anuglar momentum     */


/****************************/
/*      PHYSICAL CONSTANT   */
/****************************/

const double HBAR       =  6.58212196e-22 ; /* Planck's constant/2pi [MeV sec]*/
const double HBARSQ     =  4.33243296e-43 ; /* HBAR*HBAR                      */
const double COULOMB    =  1.60217733e-19 ; /* J = 1eV                        */
const double COULOMBSQ  =  2.56697220e-38 ; /* COULOMB*COULOMB                */
const double PERMITTIV  =  5.60958617e+37 ; /* permittivity [MeV fm /C^2]     */

const double AMUNIT     =  931.4943335    ; /* MeV = 1amu                     */
const double VLIGHT     =  2.99792458e+23 ; /* light velocty [fm/sec]         */
const double VLIGHTSQ   =  8.98755179e+46 ; /* VLIGHT*VLIGHT                  */

// const double MNEUTRON   =  1.008664891    ; /* Neutron Weight  [amu]          */
// const double MPROTON    =  1.007276487    ; /* Proton Weight   [amu]          */
// const double MDEUTERON  =  2.0141018      ; /* Deuteron Weight [amu]          */
// const double MTRITON    =  3.0160494      ; /* Triton Weight   [amu]          */
// const double MHELIUM3   =  3.0260294      ; /* Helium-3 Weight [amu]          */
// const double MALPHA     =  4.0026032      ; /* Alpha Weight    [amu]          */

const double ENEUTRON   =   8.071431      ; /* Neutron Mass Excess  [MeV]     */
const double EPROTON    =   7.289034      ; /* Proton Mass Excess   [MeV]     */
const double EDEUTERON  =  13.13584       ; /* Deuteron Mass Excess [MeV]     */
const double ETRITON    =  14.94994       ; /* Triron Mass Excess   [MeV]     */
const double EHELIUM3   =  14.93132       ; /* Helium-3 Mass Excess [MeV]     */
const double EALPHA     =   2.42492       ; /* Alpha Mass Excess    [MeV]     */
const double EELECTRON  =  0.510998910    ; /* Electron Mass        [MeV]     */

const double CSPO       =   2.04553       ; /* (HBAR/Mpi_meson/VLIGHT)^2      */
const double ALPHA      = 7.2973525376e-03; /* fine structure constant        */





#endif //__PHYSICS_H__
