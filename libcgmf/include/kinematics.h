/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file kinematics.h

  \brief Kinematic transformation from CM to LAB frames

*/
#ifndef __KINEMATICS_H__
#define __KINEMATICS_H__

#include "config-ff.h"

void boost (double fragmentMomentum[3], int Zf, int Af,
            int neutronMultiplicity,
            double cmNeutronEnergies   [MAX_NUMBER_NEUTRONS],
            double neutronEnergies     [MAX_NUMBER_NEUTRONS],
            double cmNEutronVelocities [MAX_NUMBER_NEUTRONS][3],
            double neutronVelocities   [MAX_NUMBER_NEUTRONS][3],
            int gammaMultiplicity,
            double cmGammaEnergies     [MAX_NUMBER_GAMMAS],
            double gammaEnergies       [MAX_NUMBER_GAMMAS],
            double cmGammaVelocities   [MAX_NUMBER_GAMMAS][3],
            double gammaVelocities     [MAX_NUMBER_GAMMAS][3]
            );

void boost_pfn (int Zc, int Ac,double En,double TKE,int Al,int Zl,int Ah,int Zh,
            int neutronMultiplicity_pre, double anisotropyCoefficient, bool preeq_flag,
            double cmNeutronEnergies [MAX_NUMBER_NEUTRONS_PRE], // pre-fission neutron energies in the CM (input) 
            double neutronEnergies   [MAX_NUMBER_NEUTRONS_PRE],  // pre-fission neutron energies in the lab
            int neutronType_pre      [MAX_NUMBER_NEUTRONS_PRE], // this is a placeholder, to be used later to identify
                                                            // either pre-equilibrium or statistical neutrons
            double cmNeutronVelocities [MAX_NUMBER_NEUTRONS_PRE][3],
            double neutronVelocities [MAX_NUMBER_NEUTRONS_PRE][3],
            double Pcostheta[THETA_STEPS][2],
	    double Ppreeq[THETA_STEPS][2],
            double *lightFragmentMomentum, double * heavyFragmentMomentum
            ) ;

double neutron_emission_boost(int Z_c, int A_c, double * p_c, 
			      double cmNeutronEnergy,
			      double *neutronVelocities, 
			      double *cmNeutronVelocities,
			      bool preeq_flag, double Pcostheta[THETA_STEPS][2],
			      double Ppreeq[THETA_STEPS][2]);

#endif //__KINEMATICS_H__
