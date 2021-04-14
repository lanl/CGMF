/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file kinematics.cpp

  \brief Kinematics equations for moving fragments, neutrons and gamma rays

*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <functional>

using namespace std;

#include "cgm.h"
#include "masstable.h"
#include "physics.h"
#include "kinematics.h"
#include "FissionFragments.h"

extern std::function< double(void) > rng_cgm;

void boost (double fragmentMomentum[3], int Zf, int Af,
            int    neutronMultiplicity,
            double cmNeutronEnergies [MAX_NUMBER_NEUTRONS],
            double neutronEnergies   [MAX_NUMBER_NEUTRONS],
            double cmNeutronVelocities [MAX_NUMBER_NEUTRONS][3],
            double neutronVelocities [MAX_NUMBER_NEUTRONS][3],
            int    gammaMultiplicity,
            double cmGammaEnergies   [MAX_NUMBER_GAMMAS],
            double gammaEnergies     [MAX_NUMBER_GAMMAS],
            double cmGammaVelocities [MAX_NUMBER_GAMMAS][3],
            double gammaVelocities   [MAX_NUMBER_GAMMAS][3]
            ) {
  
  double phi, cmTheta, cosCmTheta; // direction of neutron emission in c.m. of fragment
  double Mf, Mn; // fragment and neutron masses
  double Ecm; // center-of-mass neutron energy
  
  double pn; // neutron momentum in center-of-mass of fragment
  double pnx2, pny2, pnz2; // neutron momentum components (x,y,z) in c.m. frame
  double vx2, vy2, vz2; // neutron velocity components (x,y,z) in c.m. frame
  double vx, vy, vz;    // neutron velocity components (x,y,z) in lab. frame
  double v;
  
  double KEf; // fission fragment kinetic energy (MeV)
  
  // pre-neutron emission fission fragment momentum (x,y,z) in lab. frame
  double pfx = fragmentMomentum[0];
  double pfy = fragmentMomentum[1];
  double pfz = fragmentMomentum[2];
  
  double pf = sqrt(pfx*pfx+pfy*pfy+pfz*pfz);
  
  double Eg; // gamma-ray energy
  
  double pgx2, pgy2, pgz2; // photon momentum in c.m. of the fragment
  double pgx, pgy, pgz;    // photon momentum in lab
  double b, g, Egl;
  
  Mf = Af*amuMeV+getMassExcess(Zf*1000+Af);
  KEf = pf*pf/(2.0*Mf);
	
  // cout << Mf << " " << KEf << endl;
  // exit (0);

  // fragment direction
  double costhetaF = pfz/pf;
  double sinthetaF = sqrt(1.0-costhetaF*costhetaF); // abs(pfx/pf);
  double phiF;
	
  double epsilon=1e-16;
  if (pfx<0) {
    phiF = PI-atan(-pfy/pfx);
  } else if (pfy<0) {
    phiF = PI2+atan(pfy/(pfx+epsilon));
  } else {
    phiF = atan(pfy/(pfx+epsilon));
  }	

  // fragment direction, in reference frame of fragment (pfy2 and pfz2 should be zero!)
  double pfx2 = sinthetaF*cos(phiF)*pfx+sinthetaF*sin(phiF)*pfy+costhetaF*pfz; // along e_r
  double pfy2 = costhetaF*cos(phiF)*pfx + sin(phiF)*costhetaF*pfy - sinthetaF*pfz; // along e_theta
  double pfz2 = -sin(phiF)*pfx+cos(phiF)*pfy; // along e_phi
	
  Mn = NEUTRONMASS*amuMeV;

  /* IS: replace the loop over emitted neutrons with the routine neutron_emission_boost()
  
  // loop over emitted neutrons ------------------------------------------------
  for (int i=0; i<neutronMultiplicity; i++) {
    
    Mf = Af*amuMeV+mass_excess(Zf,Af);
    Ecm = cmNeutronEnergies[i];
		
    // direction of neutron emission chosen isotropic in c.m. of fission fragment
    phi = twopi*rng_cgm();
    cosCmTheta = 2.0*rng_cgm()-1.0;
    cmTheta = acos(cosCmTheta);
    
    // momentum neutron in c.m.
    pn = sqrt (2.0*Ecm*Mn);
    
    // neutron in c.m. of the fragment
    pnx2 = pn*cosCmTheta;
    pny2 = pn*sin(cmTheta)*cos(phi);
    pnz2 = pn*sin(cmTheta)*sin(phi);
    
    cmNeutronVelocities[i][0] = pnx2/Mn;
    cmNeutronVelocities[i][1] = pny2/Mn;
    cmNeutronVelocities[i][2] = pnz2/Mn;
    
    v = sqrt (pnx2*pnx2+pny2*pny2+pnz2*pnz2)/Mn;
		
    // vector addition in c.m. of fragment --> lab.
    vx2 = pnx2/Mn + pfx2/Mf;
    vy2 = pny2/Mn + pfy2/Mf;
    vz2 = pnz2/Mn + pfz2/Mf;
    
    // revert neutron vector to lab. frame
    vx = sinthetaF*cos(phiF)*vx2+costhetaF*cos(phiF)*vy2-sin(phiF)*vz2;
    vy = sinthetaF*sin(phiF)*vx2+costhetaF*sin(phiF)*vy2+cos(phiF)*vz2;
    vz = costhetaF*vx2-sinthetaF*vy2;
    
    v = sqrt (vx*vx+vy*vy+vz*vz);
    
    // save neutron velocity components
    neutronVelocities[i][0] = vx;
    neutronVelocities[i][1] = vy;
    neutronVelocities[i][2] = vz;
    
    // replace c.m. neutron energy to lab. neutron energy
    neutronEnergies[i] = 0.5*Mn*v*v;
    
    // recoil of the fragment --------------------------------------------------
    
    KEf = (Mf-Mn)/Mf*KEf + Mn/(Mf-Mn)*Ecm-2.0*sqrt(Ecm*KEf/Mf*Mn)*cosCmTheta;
    
    // new fragment velocity components
    // (pf will be multiplied later by new Mf to get back to momentum)
    pfx2 = pfx2/Mf - pnx2/(Mf-Mn);
    pfy2 = pfy2/Mf - pny2/(Mf-Mn);
    pfz2 = pfz2/Mf - pnz2/(Mf-Mn);
    
    Af--; // remove 1 neutron from evaporating fission fragment
    Mf = Af*amuMeV+mass_excess(Zf,Af);
    
    pfx2 = pfx2*Mf;
    pfy2 = pfy2*Mf;
    pfz2 = pfz2*Mf;
    
    pf = sqrt (pfx2*pfx2+pfy2*pfy2+pfz2*pfz2);
    
  } // end loop over emitted neutrons ------------------------------------------

  // post-neutron emission fragment momentum ===================================
  
  // recover pf in original (x,y,z) lab frame
  pfx = sinthetaF*cos(phiF)*pfx2+costhetaF*cos(phiF)*pfy2-sin(phiF)*pfz2;
  pfy = sinthetaF*sin(phiF)*pfx2+costhetaF*sin(phiF)*pfy2+cos(phiF)*pfz2;
  pfz = costhetaF*pfx2-sinthetaF*pfy2;
  
  pf = sqrt(pfx*pfx+pfy*pfy+pfz*pfz);
  
  // record post-neutron emission fragment momentum in the lab frame
  fragmentMomentum[0] = pfx;
  fragmentMomentum[1] = pfy;
  fragmentMomentum[2] = pfz;
  
  // ===========================================================================
  
  */

  int a_c=Af;
  for (int n=0; n<neutronMultiplicity; n++) {
      neutronEnergies[n]=neutron_emission_boost(Zf,a_c,fragmentMomentum,
						cmNeutronEnergies[n],
						neutronVelocities[n], 
						cmNeutronVelocities[n],
						false,Pcostheta,Ppreeq);
      
      a_c--; // neutron emission
    }	
  pf = 0.;
  for(int n=0;n<3;n++)
    pf+=fragmentMomentum[n]*fragmentMomentum[n];
  pf=sqrt(pf);
  // loop over gamma emissions -------------------------------------------------
  Mf = a_c*amuMeV+mass_excess(Zf,a_c);
  
  for (int i=0; i<gammaMultiplicity; i++) {
    
    // direction of gamma emission chosen isotropic in c.m. of fission fragment
    phi = twopi*rng_cgm();
    cosCmTheta = 2.0*rng_cgm()-1.0;
    cmTheta = acos(cosCmTheta);
    
    Eg = cmGammaEnergies[i];
    
    // photon momentum in c.m. of the fragment
    pgx2 = Eg*cosCmTheta;
    pgy2 = Eg*sin(cmTheta)*cos(phi);
    pgz2 = Eg*sin(cmTheta)*sin(phi);
    
    cmGammaVelocities[i][0] = pgx2;
    cmGammaVelocities[i][1] = pgy2;
    cmGammaVelocities[i][2] = pgz2;
        
    //		if (index!=0) cmTheta=PI-cmTheta;
    b = pf/Mf; // why is this not outside?
    Egl= Eg*sqrt(1.0-b*b)/(1.0-b*cos(cmTheta)); // Doppler shift
    
    // Boost photon in c.m. of fragment --> lab
    
    // FOR NOW, keep same orientation but change gamma-ray energy only
    pgx = Egl*cosCmTheta;
    pgy = Egl*sin(cmTheta)*cos(phi);
    pgz = Egl*sin(cmTheta)*sin(phi);
    
    // save gamma velocity components
    gammaVelocities[i][0] = pgx;
    gammaVelocities[i][1] = pgy;
    gammaVelocities[i][2] = pgz;
    
    // replace c.m. energy with Doppler shifted energy
    gammaEnergies[i] = Egl;
    
  } // end loop over gamma emissions -------------------------------------------
  
  
  return;
  
}

/*==========================================================================
This function takes as input the initial momentum of a nucleus,
and the neutron CM energy, and computes the kinematics, 
returning the lab neutron velocity and the momentum of the daughter nucleus
Maybe we should overload this function, so that we don't check for
pre-equlibrium neutron if we know it's statistical
p_c contains the post neutron emission momentum, if the initial
one is desired, one needs to save it
 */
double neutron_emission_boost(int Z_c, int A_c, double * p_c, 
			      double cmNeutronEnergy,
			      double *neutronVelocities, 
			      double *cmNeutronVelocities,
			      bool preeq_flag, double Pcostheta[THETA_STEPS][2],
			      double Ppreeq[THETA_STEPS][2]){

  double rnum,x1,y1,x2,y2;
  double cosTheta, sinTheta;

  double M_c = A_c*amuMeV+getMassExcess(Z_c*1000+A_c);  
  double Mn = NEUTRONMASS*amuMeV;

  double phi = twopi*rng_cgm();
  if (preeq_flag){
    // sample from the correct angular distribution
    rnum = rng_cgm();
    x1 = 0;
    x2 = 0;
    y1 = 0;
    y2 = 0;
    for (int j=0;j<THETA_STEPS-1;j++){
      if (rnum >= Ppreeq[j][1]){
        x1 = Ppreeq[j][0];
        y1 = Ppreeq[j][1];
        x2 = Ppreeq[j+1][0];
        y2 = Ppreeq[j+1][1];
      } else { break; }
    }
    cosTheta = ((x2-x1)*rnum - y1*x2 + y2*x1)/(y2-y1); // interpolate between tabulated CDF points
  } else {
    cosTheta = 2.0*rng_cgm()-1.0;
  }
  sinTheta = sqrt(1.-cosTheta*cosTheta);

  double vn_cm = sqrt(2.*cmNeutronEnergy/Mn);
  cmNeutronVelocities[0]=vn_cm*sinTheta*cos(phi);
  cmNeutronVelocities[1]=vn_cm*sinTheta*sin(phi);
  cmNeutronVelocities[2]=vn_cm*cosTheta;
  double neutronEnergy=0.;
  for(int i=0;i<3;i++){
    neutronVelocities[i]=cmNeutronVelocities[i]+p_c[i]/M_c; // lab velocity
    neutronEnergy+=neutronVelocities[i]*neutronVelocities[i];
    p_c[i] -= Mn*neutronVelocities[i];
  }
  return neutronEnergy*.5*Mn;
}


/*===========================================================================
Computes the boost of the pre-fission neutrons and the momenta of the two
fragments after fission, based on non-relativistic kinematics
Note that Z_c, A_c, En characterize the initial compound system
=============================================================================*/
void boost_pfn (int Z_c, int A_c,double En,double TKE,int A_l,int Z_l,int A_h,int Z_h,
    int neutronMultiplicity_pre, double anisotropyCoefficient, bool preeq_flag,
    double cmNeutronEnergies [MAX_NUMBER_NEUTRONS_PRE], // pre-fission neutron energies in the CM (input) 
    double neutronEnergies   [MAX_NUMBER_NEUTRONS_PRE],  // pre-fission neutron energies in the lab
    int neutronType_pre      [MAX_NUMBER_NEUTRONS_PRE], // this is a placeholder, to be used later to identify
                                                    // either pre-equilibrium or statistical neutrons
    double cmNeutronVelocities [MAX_NUMBER_NEUTRONS_PRE][3],
    double neutronVelocities [MAX_NUMBER_NEUTRONS_PRE][3],
    double Pcostheta[THETA_STEPS][2],
    double Ppreeq[THETA_STEPS][2],
    double *lightFragmentMomentum, double *heavyFragmentMomentum
            ) {
  double M_c, Mn; // compound nucleus mass and neutron masses
  double phi, cosTheta, sinTheta; // direction of neutron emission in c.m. of fragment
  double Ecm; // center-of-mass neutron energy
  double rnum; // random number, saved for comparison purpose
  double x1,x2,y1,y2; // interpolate between the tabulated values of Pcostheta

  double p_c[3];
  int a_c=A_c;  // running compound system mass taking into account the neutron emission

  Mn = NEUTRONMASS*amuMeV;
  p_c[0]=0.;
  p_c[1]=0.;
  p_c[2]=sqrt(2.*Mn*En); // maybe we should use the reduced mass??
  int n0=0;
  if(preeq_flag){
    neutronEnergies[0]=neutron_emission_boost(Z_c,a_c,p_c,
					      cmNeutronEnergies[0],
					      neutronVelocities[0], 
					      cmNeutronVelocities[0],
					      true,Pcostheta,Ppreeq);
    n0=1;
  }      
    
  for (int n=n0;n<neutronMultiplicity_pre;n++){ 
      neutronEnergies[n]=neutron_emission_boost(Z_c,a_c,p_c,
						cmNeutronEnergies[n],
						neutronVelocities[n], 
						cmNeutronVelocities[n],
						false,Pcostheta,Ppreeq);
      
      a_c--; // neutron emission
  }
  M_c = a_c*amuMeV+getMassExcess(Z_c*1000+a_c);  

  // All the prefission neutrons emitted, now the light fragment is emitted isotropically 
  // (unless anisotropy coefficient available in anisotropy.dat)
  // and the heavy-fragment momentum is calculated from the momentum conservation

  // first, the TKE in the CM is given by TKE-.5*(M_L+M_H)*v_c^2

  double M_l=A_l*amuMeV+getMassExcess(Z_l*1000+A_l);
  double M_h=A_h*amuMeV+getMassExcess(Z_h*1000+A_h);

  double v_c2=0.;
  /* IS: This is be the correction due to the boost */
  for(int i=0;i<3;i++)
    v_c2+=p_c[i]*p_c[i];
  v_c2/=(M_c*M_c);

  double TKE_cm=TKE-.5*(M_l+M_h)*v_c2; // could this be negative?
  double Pf0 = sqrt(2.*TKE_cm/(1./M_l+1/M_h)); // magnitude of the momentum for each fragment
  phi = twopi*rng_cgm();

  // AL: Including anisotropy in the emission of the fission fragments, for select nuclei
  if (anisotropyCoefficient==1.0){
    cosTheta = 2.0*rng_cgm()-1.0; // sample isotropic when 
  }
  else{
    // sample from the cos(theta) CDF distribution calculated in computeFragmentAngularDistribution
    rnum = rng_cgm();
    for (int j=0;j<THETA_STEPS-1;j++){
      if (rnum >= Pcostheta[j][1]){
        x1 = Pcostheta[j][0];
        y1 = Pcostheta[j][1];
        x2 = Pcostheta[j+1][0];
        y2 = Pcostheta[j+1][1];
      } else { break; }
    }
    cosTheta = ((x2-x1)*rnum - y1*x2 + y2*x1)/(y2-y1); // interpolate between tabulated CDF points
  }
  sinTheta = sqrt(1.0 - cosTheta*cosTheta);

  lightFragmentMomentum[0]=Pf0*cos(phi)*sinTheta+M_l*p_c[0]/M_c;
  lightFragmentMomentum[1]=Pf0*sin(phi)*sinTheta+M_l*p_c[1]/M_c;
  lightFragmentMomentum[2]=Pf0*cosTheta         +M_l*p_c[2]/M_c;
  
  for(int i=0;i<3;i++){
    heavyFragmentMomentum[i]=p_c[i]-lightFragmentMomentum[i];
  }
    
  return;
  
}
