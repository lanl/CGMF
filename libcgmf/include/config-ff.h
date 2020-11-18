/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file config-ff.h

  \brief Additional configuration file for CGMF

*/


#ifndef __CONFIG_FF_H__
#define __CONFIG_FF_H__

#ifdef __cplusplus
#include <string>
#endif

// MER: Should this be a command line option or something?
#define HISTORIES // to write out Monte Carlo history file
//#undef HISTORIES

/*----------------------------------------------------------
 Prompt gamma-ray results are saved into a format usable for
 GEANT simulations of the DANCE g-ray detector.
 */

//#define GEANT
#undef GEANT

//const std::string WORKDIR = "./";  -->> Moved to FissionFragments.h

const int    NUMA   =  300; // number of masses A
const int    NUMZ   =  100; // number of charges Z
const int    NUMTKE =  300; // number of Total Kinetic Energy values

const int      NUME =  401; // number of energies in level density tables; dE=0.25 MeV; up to Emax=100 MeV
const double deltaE = 0.25; // energy-bin size used in level density tables
const int     NUME2 = 32;

const int    NUMdZ  =   21; // [-dZ:+dZ] if dZ=10 for charge distribution around most probable Zp[A] 

const int   NUMMULT =   50; // number of multiplicities

const int NUMANGLES =   73; // number of angles in angular distribution
const double dTheta =  2.5; // angular bins (degrees)

const int MAX_NUMBER_PARTICLES = 50; // max. number of particles (n or g) emitted per fragment in a fission event

#define   MAX_NUMBER_NEUTRONS 10
#define   MAX_NUMBER_GAMMAS   50
const int MAX_NUMBER_NEUTRONS_PRE = 3; // max. number of pre-fission neutrons emitted

const int NUMBER_SPECTRUM_ENERGY_GRID = 641; //551;

const int THETA_STEPS = 50; // number of steps in the cos(theta) CDF distribution

// to pass as argument to MCNP fission event
struct fissionEventType {

	int nu = 0;  // number of neutrons
	double neutronEnergies[MAX_NUMBER_NEUTRONS] = {};
	double neutronDircosu[MAX_NUMBER_NEUTRONS] = {};
	double neutronDircosv[MAX_NUMBER_NEUTRONS] = {};
	double neutronDircosw[MAX_NUMBER_NEUTRONS] = {};
	double neutronAges[MAX_NUMBER_NEUTRONS] = {};
	
	double cmNeutronEnergies[MAX_NUMBER_NEUTRONS] = {};
	double cmNeutronDircosu[MAX_NUMBER_NEUTRONS] = {};
	double cmNeutronDircosv[MAX_NUMBER_NEUTRONS] = {};
	double cmNeutronDircosw[MAX_NUMBER_NEUTRONS] = {};
	
	int nug= 0; // number of photons
	double photonEnergies[MAX_NUMBER_GAMMAS] = {};
	double photonDircosu[MAX_NUMBER_GAMMAS] = {};
	double photonDircosv[MAX_NUMBER_GAMMAS] = {};
	double photonDircosw[MAX_NUMBER_GAMMAS] = {};
	double photonAges[MAX_NUMBER_GAMMAS] = {};
	
	double cmPhotonEnergies[MAX_NUMBER_GAMMAS] = {};
	double cmPhotonDircosu[MAX_NUMBER_GAMMAS] = {};
	double cmPhotonDircosv[MAX_NUMBER_GAMMAS] = {};
	double cmPhotonDircosw[MAX_NUMBER_GAMMAS] = {};

	int npfn = 0; // number of pre-fission neutrons
	
	double preFragmentMomentum[3] = {}; // before neutron emission
	double postFragmentMomentum[3] = {}; // after neutron emission

	double fragmentDircosu = 0.0;
        double fragmentDircosv = 0.0;
        double fragmentDircosw = 0.0;
	
	int A = 0, Z = 0; // fission fragment mass and charge numbers
	double KE = 0.0, KEpost = 0.0; // fragment kinetic energy (pre- and post-neutron emission) 
	double U = 0.0; // initial excitation energy
	float spin = 0.0; // initial spin
	int parity = 0; // initial parity

        fissionEventType() = default;
	
};

#endif //__CONFIG_FF_H__
