/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf_help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file cgmfEvents_capi.cpp

  \brief CGM and CGMF C-API
 
 */

#include <string>

#include "cgmfEvents_capi.h"
#include "cgmfEvents.h"

using namespace std;

cgmEvent*  cgm_event;
cgmfEvent* cgmf_fevent;
cgmfYields* cgmf_yields;

/*!
\brief Set the random number generator used throughout the code
*/  
void setrngdptr(double (*funcptr) (void)) {
  set_rng(funcptr);
}

/*!
\brief Used to check the directory path containing auxiliary data files, usually in the main driver file.
*/
int checkdatapath (const char* datapath) {
  return checkdatapath(string(datapath));
}

/*!
\brief Used to set the directory path containing auxiliary data files, usually in the main driver file.
*/
void setdatapath (const char* datapath) {
  setdatapath(string(datapath));
}

/*!
\brief Instantiates a new cgm event.
*/
void cgm(int dza, double dex, int *dmcl_nlines, double dmcl_glines[30],
         int *dmcl_nlinesn, double dmcl_elinesn[30]) {
  if (cgm_event != 0) delete cgm_event;

  cgm_event = new cgmEvent(dza,dex,0.0);

  *dmcl_nlines = (*cgm_event).getPhotonNu();
  for (int ii=0; ii<*dmcl_nlines; ii++) dmcl_glines[ii] = (*cgm_event).getPhotonEnergy(ii);
  *dmcl_nlinesn = (*cgm_event).getNeutronNu();
  for (int ii=0; ii<*dmcl_nlinesn; ii++) dmcl_elinesn[ii] = (*cgm_event).getNeutronEnergy(ii);
}

/*!
\brief Instantiates a new element of the cgmfEvent class.
Given a charge (Z), mass (A), and incident particle energy (erg), it instantiates a new cgmfEvent.
*/
void cgmf_genfissevent(int za, double erg, double tme, double tmw) {
  if (cgmf_fevent != 0) delete cgmf_fevent;
  cgmf_fevent = new cgmfEvent(za,erg,tme,tmw);
}

/*!
\brief Generates pre-neutron emission fission yields in mass, charge, kinetic energy, excitation energy,
spin and parity. 

\param za: ZAID of target nucleus
\param erg: energy of the incident particle in the laboratory frame, in MeV
\param nevents: number of fission events to generate

*/
void cgmf_genfissyields (int za, double erg, int nevents, string outfilename) {
  cgmf_yields = new cgmfYields (za, erg, nevents, outfilename);
}

/*!
\brief Returns the neutron multiplicity of the fission event
*/	
int cgmf_getnnu() {
  return cgmf_fevent->getNeutronNu();
}

/*!
\brief Returns the lab and c.o.m. energies of the i^th emitted neutron
*/
double cgmf_getnerg(int i) { return cgmf_fevent->getNeutronEnergy(i); }
double cgmf_getcmnerg(int i) { return cgmf_fevent->getCmNeutronEnergy(i); }

/*!
\brief Returns the first (u) directional cosine in the laboratory frame of the i^th emitted neutron
*/
double cgmf_getndircosu(int i) { return cgmf_fevent->getNeutronDircosu(i); }
double cgmf_getcmndircosu(int i) { return cgmf_fevent->getCmNeutronDircosu(i); }

/*!
\brief Returns the second (v) directional cosine in the laboratory frame of the i^th emitted neutron
*/
double cgmf_getndircosv(int i) { return cgmf_fevent->getNeutronDircosv(i); }
double cgmf_getcmndircosv(int i) { return cgmf_fevent->getCmNeutronDircosv(i); }

/*!
\brief Returns the third (w) directional cosine in the laboratory frame of the i^th emitted neutron
*/
double cgmf_getndircosw(int i) { return cgmf_fevent->getNeutronDircosw(i); }
double cgmf_getcmndircosw(int i) { return cgmf_fevent->getCmNeutronDircosw(i); }

/*!
\brief Returns the time of emission (age) since fission of the i^th emitted neutron
*/
double cgmf_getntme(int i) {
  return cgmf_fevent->getNeutronAge(i);
}

/*!
\brief Returns the gamma multiplicity of the fission event
*/
int cgmf_getgnu() {
  return cgmf_fevent->getPhotonNu();
}

/*!
\brief Returns the lab energy of the i^th emitted photon
*/
double cgmf_getgerg(int i) {
  return cgmf_fevent->getPhotonEnergy(i);
}

/*!
\brief Returns the first (u) directional cosine in the laboratory frame of the i^th emitted photon
*/
double cgmf_getgdircosu(int i) {
  return cgmf_fevent->getPhotonDircosu(i);
}

/*!
\brief Returns the second (v) directional cosine in the laboratory frame of the i^th emitted photon
*/
double cgmf_getgdircosv(int i) {
  return cgmf_fevent->getPhotonDircosv(i);
}

/*!
\brief Returns the third (w) directional cosine in the laboratory frame of the i^th emitted photon
*/
double cgmf_getgdircosw(int i) {
  return cgmf_fevent->getPhotonDircosw(i);
}

/*!
\brief Returns the time of emission (age) since fission of the i^th emitted photon
*/
double cgmf_getgtme(int i) {
  return cgmf_fevent->getPhotonAge(i);
}

/*
\brief Delete fission events and clean up arrays
*/
void cgmf_clean() {
  if (cgmf_fevent != 0) delete cgmf_fevent;
//  riplDiscreteLevelsCleanup();
}

int cgmf_getlfmass() { return cgmf_fevent->getLightFragmentMass(); }
int cgmf_gethfmass() { return cgmf_fevent->getHeavyFragmentMass(); }

int cgmf_getlfcharge() { return cgmf_fevent->getLightFragmentCharge(); }
int cgmf_gethfcharge() { return cgmf_fevent->getHeavyFragmentCharge(); }

double cgmf_getlfke() { return cgmf_fevent->getLightFragmentKineticEnergy(); }
double cgmf_gethfke() { return cgmf_fevent->getHeavyFragmentKineticEnergy(); }

double cgmf_getlfkepost() { return cgmf_fevent->getLightFragmentKineticEnergyPost(); }
double cgmf_gethfkepost() { return cgmf_fevent->getHeavyFragmentKineticEnergyPost(); }

double cgmf_getlfu() { return cgmf_fevent->getLightFragmentExcitationEnergy(); }
double cgmf_gethfu() { return cgmf_fevent->getHeavyFragmentExcitationEnergy(); }

float cgmf_getlfspin() { return cgmf_fevent->getLightFragmentSpin(); }
float cgmf_gethfspin() { return cgmf_fevent->getHeavyFragmentSpin(); }

int cgmf_getlfparity() { return cgmf_fevent->getLightFragmentParity(); }
int cgmf_gethfparity() { return cgmf_fevent->getHeavyFragmentParity(); }

int cgmf_getlfnnu() { return cgmf_fevent->getLightFragmentNeutronNu(); }
int cgmf_gethfnnu() { return cgmf_fevent->getHeavyFragmentNeutronNu(); }

// pre-fission neutron number and directional cosines
int    cgmf_getprennu()           { return cgmf_fevent->getPreFissionNeutronNu(); }
double cgmf_getprenerg(int i)     { return cgmf_fevent->getPreFissionNeutronEnergy(i); }
double cgmf_getprendircosu(int i) { return cgmf_fevent->getPreFissionNeutronDircosu(i); }
double cgmf_getprendircosv(int i) { return cgmf_fevent->getPreFissionNeutronDircosv(i); }
double cgmf_getprendircosw(int i) { return cgmf_fevent->getPreFissionNeutronDircosw(i); }


int cgmf_getlfgnu() { return cgmf_fevent->getLightFragmentPhotonNu(); }
int cgmf_gethfgnu() { return cgmf_fevent->getHeavyFragmentPhotonNu(); }

double cgmf_getlfdircosu() { return cgmf_fevent->getLightFragmentDircosu(); }
double cgmf_gethfdircosu() { return cgmf_fevent->getHeavyFragmentDircosu(); }

double cgmf_getlfdircosv() { return cgmf_fevent->getLightFragmentDircosv(); }
double cgmf_gethfdircosv() { return cgmf_fevent->getHeavyFragmentDircosv(); }

double cgmf_getlfdircosw() { return cgmf_fevent->getLightFragmentDircosw(); }
double cgmf_gethfdircosw() { return cgmf_fevent->getHeavyFragmentDircosw(); }

double cgmf_getlfpremomentum_x () { return cgmf_fevent->getLightFragmentPreMomentumX (); }
double cgmf_gethfpremomentum_x () { return cgmf_fevent->getHeavyFragmentPreMomentumX (); }

double cgmf_getlfpremomentum_y () { return cgmf_fevent->getLightFragmentPreMomentumY (); }
double cgmf_gethfpremomentum_y () { return cgmf_fevent->getHeavyFragmentPreMomentumY (); }

double cgmf_getlfpremomentum_z () { return cgmf_fevent->getLightFragmentPreMomentumZ (); }
double cgmf_gethfpremomentum_z () { return cgmf_fevent->getHeavyFragmentPreMomentumZ (); }

double cgmf_getlfpostmomentum_x () { return cgmf_fevent->getLightFragmentPostMomentumX (); }
double cgmf_gethfpostmomentum_x () { return cgmf_fevent->getHeavyFragmentPostMomentumX (); }

double cgmf_getlfpostmomentum_y () { return cgmf_fevent->getLightFragmentPostMomentumY (); }
double cgmf_gethfpostmomentum_y () { return cgmf_fevent->getHeavyFragmentPostMomentumY (); }

double cgmf_getlfpostmomentum_z () { return cgmf_fevent->getLightFragmentPostMomentumZ (); }
double cgmf_gethfpostmomentum_z () { return cgmf_fevent->getHeavyFragmentPostMomentumZ (); }

