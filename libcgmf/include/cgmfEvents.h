/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file cgmfEvents.h

  \brief CGM and CGMF classes

*/


#ifndef __CGMFEVENTS_H__
#define __CGMFEVENTS_H__

#ifdef __cplusplus
#include <string>
#endif
#include <vector>
#include <functional>

#include "config-ff.h"

int  checkdatapath(std::string);
void setdatapath(std::string);
void cgmf_cleanup(void);

void set_rng(std::function<double(void)>);
void set_rng(double (*) (void));


/* cgmEvent class: base class for cgmfEvent */

class cgmEvent {

public:

	cgmEvent ();
	cgmEvent (int , double , double);
	~cgmEvent ();
	
	// neutrons ------------------------------------------------------------------
	int getNeutronNu () { return neutronNu; }
	
	double getNeutronEnergy(int index) {
		if (index >= 0 && index < neutronNu) return neutronEnergies[index];
		else return -1;
	}
	double getNeutronDircosu(int index) {
		if (index >= 0 && index < neutronNu) return neutronDircosu[index];
		else return -1;
	}
	double getNeutronDircosv(int index) {
		if (index >= 0 && index < neutronNu) return neutronDircosv[index];
		else return -1;
	}
	double getNeutronDircosw(int index) {
		if (index >= 0 && index < neutronNu) return neutronDircosw[index];
		else return -1;
	}
	double getNeutronAge(int index) {
		if (index >= 0 && index < neutronNu) return neutronAges[index];
		else return -1;
	}

	// neutrons (in center-of-mass of emitting fragment) -------------------------

	double getCmNeutronEnergy(int index) {
		if (index >= 0 && index < neutronNu) return cmNeutronEnergies[index];
		else return -1;
	}
	double getCmNeutronDircosu(int index) {
		if (index >= 0 && index < neutronNu) return cmNeutronDircosu[index];
		else return -1;
	}
	double getCmNeutronDircosv(int index) {
		if (index >= 0 && index < neutronNu) return cmNeutronDircosv[index];
		else return -1;
	}
	double getCmNeutronDircosw(int index) {
		if (index >= 0 && index < neutronNu) return cmNeutronDircosw[index];
		else return -1;
	}
	
	// photons -------------------------------------------------------------------
	int getPhotonNu  () { return photonNu; }
	
	double getPhotonEnergy(int index) {
		if (index >= 0 && index < photonNu) return photonEnergies[index];
		else return -1;
	}
	double getPhotonDircosu(int index) {
		if (index >= 0 && index < photonNu) return photonDircosu[index];
		else return -1;
	}
	double getPhotonDircosv(int index) {
		if (index >= 0 && index < photonNu) return photonDircosv[index];
		else return -1;
	}
	double getPhotonDircosw(int index) {
		if (index >= 0 && index < photonNu) return photonDircosw[index];
		else return -1;
	}
	double getPhotonAge(int index) {
		if (index >= 0 && index < photonNu) return photonAges[index];
		else return -1;
	}

protected:

	int neutronNu;
	double* neutronEnergies;
	double* neutronDircosu;
	double* neutronDircosv;
	double* neutronDircosw;
	double* neutronAges;

	double* cmNeutronEnergies;
	double* cmNeutronDircosu;
	double* cmNeutronDircosv;
	double* cmNeutronDircosw;
	
	int photonNu;
	double* photonEnergies;
	double* photonDircosu;
	double* photonDircosv;
	double* photonDircosw;
	double* photonAges;

	void outputOptions (unsigned int);
	void allocateMemory ();
    int cgmCheckRange(int, int, double);

private:
	
	void initialization (unsigned int);

};

/* cgmfEvent class: derived from cgmEvent base class */

class cgmfEvent: public cgmEvent {

public:
	
	cgmfEvent (int , double , double, double = -1.0);
	cgmfEvent ();
	~cgmfEvent ();

	// Reaction -----------------------------------------------------------------
	int    getTargetNucleus () { return ZAIDt; }
	double getIncidentEnergy () { return incidentEnergy; }

	// Fragments ----------------------------------------------------------------
	int    getLightFragmentMass ()   { return eventLF.A; }
	int    getLightFragmentCharge () { return eventLF.Z; }
	double getLightFragmentKineticEnergy () { return eventLF.KE; }
	double getLightFragmentKineticEnergyPost () { return eventLF.KEpost; }
	double getLightFragmentExcitationEnergy () { return eventLF.U; }
	float  getLightFragmentSpin () { return eventLF.spin; }
	int    getLightFragmentParity () { return eventLF.parity; }
	int    getLightFragmentNeutronNu () { return eventLF.nu; }
	int    getLightFragmentPhotonNu () { return eventLF.nug; }

	double getLightFragmentPreMomentumX () { return eventLF.preFragmentMomentum[0]; }
	double getLightFragmentPreMomentumY () { return eventLF.preFragmentMomentum[1]; }
	double getLightFragmentPreMomentumZ () { return eventLF.preFragmentMomentum[2]; }

	double getLightFragmentDircosu () { return eventLF.fragmentDircosu; }
	double getLightFragmentDircosv () { return eventLF.fragmentDircosv; }
	double getLightFragmentDircosw () { return eventLF.fragmentDircosw; }

	double getLightFragmentPostMomentumX () { return eventLF.postFragmentMomentum[0]; }
	double getLightFragmentPostMomentumY () { return eventLF.postFragmentMomentum[1]; }
	double getLightFragmentPostMomentumZ () { return eventLF.postFragmentMomentum[2]; }

	int    getHeavyFragmentMass ()   { return eventHF.A; }
	int    getHeavyFragmentCharge () { return eventHF.Z; }
	double getHeavyFragmentKineticEnergy () { return eventHF.KE; }
	double getHeavyFragmentKineticEnergyPost () { return eventHF.KEpost; }
	double getHeavyFragmentExcitationEnergy () { return eventHF.U; }
	float  getHeavyFragmentSpin () { return eventHF.spin; }
	int    getHeavyFragmentParity () { return eventHF.parity; }
	int    getHeavyFragmentNeutronNu () { return eventHF.nu; }
	int    getHeavyFragmentPhotonNu () { return eventHF.nug; }

	double getHeavyFragmentPreMomentumX () { return eventHF.preFragmentMomentum[0]; }
	double getHeavyFragmentPreMomentumY () { return eventHF.preFragmentMomentum[1]; }
	double getHeavyFragmentPreMomentumZ () { return eventHF.preFragmentMomentum[2]; }

	double getHeavyFragmentDircosu () { return eventHF.fragmentDircosu; }
	double getHeavyFragmentDircosv () { return eventHF.fragmentDircosv; }
	double getHeavyFragmentDircosw () { return eventHF.fragmentDircosw; }

	double getHeavyFragmentPostMomentumX () { return eventHF.postFragmentMomentum[0]; }
	double getHeavyFragmentPostMomentumY () { return eventHF.postFragmentMomentum[1]; }
	double getHeavyFragmentPostMomentumZ () { return eventHF.postFragmentMomentum[2]; }

 	int    getPreFissionNeutronNu () { return preFissionNeutronNu; }

 	double getPreFissionNeutronEnergy(int index) {
 		if (index >= 0 && index < preFissionNeutronNu) return preFissionNeutronEnergies[index];
 		else return -1;
 	}

 	double getPreFissionNeutronDircosu(int index) { 
		if (index >= 0 && index < preFissionNeutronNu) return preFissionNeutronDircosu[index];
		else return -1;
	}
 	double getPreFissionNeutronDircosv(int index) { 
		if (index >= 0 && index < preFissionNeutronNu) return preFissionNeutronDircosv[index];
		else return -1;
	}
 	double getPreFissionNeutronDircosw(int index) { 
		if (index >= 0 && index < preFissionNeutronNu) return preFissionNeutronDircosw[index];
		else return -1;
	}

protected:

//	cgmfEventType lfevent, hfevent; // light and heavy fragment data
 
 	fissionEventType eventLF, eventHF;

 	int ZAIDt; // ZAID of target nucleus
 	double incidentEnergy;

	int preFissionNeutronNu;
	double* preFissionNeutronEnergies;
	double* preFissionNeutronDircosu;
	double* preFissionNeutronDircosv;
	double* preFissionNeutronDircosw;
	double* preFissionNeutronAges;

	void allocateMemory();

private:
	
	void initialization ();
	void checkInput (int, double);

};

/* cgmfYields class: derived from cgmfEvent base class */

class cgmfYields: public cgmfEvent {
	
public:
	
	cgmfYields (int, double, int, std::string);
	~cgmfYields ();
	
};

#endif //__CGMFEVENTS_H__
