/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file cgmfEvents.cpp

  \brief CGM and CGMF Classes
 
  \mainpage The CGMF Code

  CGMF implements the Hauser-Feshbach statistical theory of compound nuclear 
  reactions to the de-excitation of fission fragments by emission of prompt 
  neutrons and gamma rays.

 */

#include <sys/types.h>
#include <sys/stat.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <ctime>
//#include <netdb.h>
//#include <unistd.h>

using namespace std;

#include "cgmfEvents.h"
#include "cgm.h"
#include "global.h"
#include "terminate.h"
#include "config.h"
#include "config-ff.h"
#include "FissionFragments.h"
#include "physics.h"
#include "kcksyst.h"
#include "ripl2levels.h"
#include "rngcgm.h"
#include "cgmf_config.h"

static bool firstcalc = true;
string datadir = "";

void set_time_coincidence_window(double timew);


std::function< double(void) > rng_cgm;
void set_rng(std::function<double(void)> rng) {
  rng_cgm = rng;
}
void set_rng(double (*rng) (void)) {
  rng_cgm = rng;
}

/*!
\brief Used to check the directory path containing auxiliary data files, usually in the main driver file.
*/
int checkdatapath (string datapath) {
    struct stat info;

    if(stat( datapath.c_str(), &info ) != 0)
        return 0;
    else if(info.st_mode & S_IFDIR)
        return 1;
    else
        return 0;
}

/*!
\brief Used to set the directory path containing auxiliary data files, usually in the main driver file.
*/
void setdatapath (string path) {

  if (path == "") {
    if (checkdatapath(INSTALL_DATADIR)) {
      path = INSTALL_DATADIR;
    } else if (checkdatapath(BUILD_DATADIR)) {
      path = BUILD_DATADIR;
    } else {
      cerr << "Cannot find valid path to CGMF data ... returning" << endl;
      exit(-1);
    }
  }
  if( path.back() != '/' ) path += "/";

  datadir = path;
}


int MAX_ENERGY_BIN;   ///< Maximum number of energy bins used in the continuum
int MAX_J;            ///< Maximum number of angular momentum values in the spin distributions 
int MAX_LEVELS;       ///< Maximum number of discrete levels used in any nucleus

double ENERGY_BIN;          ///< size of the energy bins in the continuum (in MeV)
double CONTINUUM_LOWER_CUT; ///< lower-energy threshold to eliminate very small gamma energies in the continuum

int mcl_nlines;
int mcl_nlinesn;
double mcl_glines[30]={1.0};
double mcl_elinesn[30]={1.0};

Calculation   ctl;
Nucleus      *ncl;                 ///< 0: parent, 1 - MAX_COMPOUND-1: daughters etc
double       *spc[SPECTRA_OUTPUT]; ///< 0: gamma, 1: neutron, 2: electron, 3: neutrino

ostringstream oserr;


/************************************************************************************
                cgmEvent: base class constructor, destructor & methods 
 ************************************************************************************/
cgmEvent::cgmEvent () {
  return;
}


cgmEvent::cgmEvent (int isotope, double eng, double time) {

  double        exciE = 0.0; // nuclear excitation energy
  double        initJ = -1.; // initial state spin, J
  double        targE = 0.0; // target kinetic energy, when moving
  double        spinf = 0.0; // multiplication factor for initial spin distribution (sigma)
  int           initP = 0;   // initial state parity, -1 or 1
  int           targZ = 0;   // target Z number
  int           targA = 0;   // target A number
  unsigned long nsim  = 1;   // number of Monte Carlo sampling

  int          p;
  unsigned int o=0xffff,c=0;

  neutronNu = 0;
  photonNu  = 0;

  targZ = int(isotope/1000.0);
  targA = isotope-1000*targZ;
  exciE = eng;
  c = c | CALC_MC;  // hard code the -m switch

  if(firstcalc) {
    initialization(o);
    firstcalc = false;
  }
  cgmAllocateMemory();

  //---------------------------------------
  //      Check Z, A, and E Range, spin and parity
  if(c != 0){
    ctl.calc_montecarlo  = c & CALC_MC;
    ctl.calc_entrance    = c & CALC_TRAN;
  }
  if(initJ >= 0.0) ctl.init_pop =  SINGLE;
  else{
    ctl.init_pop = (ctl.calc_entrance) ? TRANSMISSION : LEVDEN;
  }

  ncl[0].za.setZA(targZ,targA);
  ncl[0].max_energy = exciE;

  if(ctl.calc_montecarlo) specMCMain(initJ,initP,targE,spinf,nsim,spc);
  else                    specMain(initJ,initP,spinf,spc,0.0);

  /*** Spectrum in the LAB frame */
  if(targE > 0.0) cgmLabSpectrum(targE/targA,ncl[0].binwidth,spc);

  // assign vars to pass back to fortran code
  cgmGetSpectra(ncl[0].de,spc);
  photonNu = mcl_nlines;
  neutronNu = mcl_nlinesn;
  allocateMemory();
  if (photonNu>0) {
    for (int ii=0; ii<mcl_nlines; ii++) photonEnergies[ii] = mcl_glines[ii];
  }
  if (neutronNu>0) {
    for (int ii=0; ii<mcl_nlinesn; ii++) neutronEnergies[ii] = mcl_elinesn[ii];
  }

  return;
}


cgmEvent::~cgmEvent () {
  if (neutronNu > 0) {
    delete [] neutronEnergies;
    delete [] neutronDircosu;
    delete [] neutronDircosv;
    delete [] neutronDircosw;
    delete [] neutronAges;
    delete [] cmNeutronEnergies;
    delete [] cmNeutronDircosu;
    delete [] cmNeutronDircosv;
    delete [] cmNeutronDircosw;
  }
  if (photonNu > 0) {
    delete [] photonEnergies;
    delete [] photonDircosu;
    delete [] photonDircosv;
    delete [] photonDircosw;
    delete [] photonAges;
  }
  cgmDeleteAllocated();
  
}
/*****************************************************************
           end of cgmEvent constructors and destructors 
 *****************************************************************/


/**********************************************************/
/*      Initializing cgmEvent                             */
/**********************************************************/
void cgmEvent::initialization (unsigned int o) {

  MAX_ENERGY_BIN      = 2500; //1000; // 2000;
  MAX_J               = 50; //50; // 70;
  MAX_LEVELS          = 200; //200; // 700;
  ENERGY_BIN          = 0.05; //0.1; // 0.05;
  CONTINUUM_LOWER_CUT = 0.02;

  outputOptions(o);

  readDiscreteLevelData();
  readkcksystdat();
  readMasses();
  readDeformations();
  readTemperatures();
  readLDP();
  readAnisotropy();
  readPreEqAngularDistributionParameters();
  readRTA();
  readSpinScaling();
  readPreeqData();

}


/**********************************************************/
/*      Setting Output Options                            */
/**********************************************************/
void cgmEvent::outputOptions(unsigned int p)
{
	ctl.print_spectrum     = DEFAULT_CAL & OUT_SPEC;
	ctl.print_init_pop     = DEFAULT_CAL & OUT_INIPOP;
	ctl.print_history      = DEFAULT_CAL & OUT_HIST;
	ctl.print_indivspec    = DEFAULT_CAL & OUT_INDIVSPEC;
	
	if(p != 0xffff){
		ctl.print_spectrum   = p & OUT_SPEC;
		ctl.print_init_pop   = p & OUT_INIPOP;
		ctl.print_history    = p & OUT_HIST;
		ctl.print_broadened  = p & OUT_BROADENED;
		ctl.print_indivspec  = p & OUT_INDIVSPEC;
	}
	
	if(ctl.print_broadened) ctl.print_spectrum = true;
}


/**********************************************************/
/*      Allocate Memory for CGM & CGMF                    */
/**********************************************************/
void cgmEvent::allocateMemory()
{
  if (neutronNu > 0) {
    neutronEnergies = new double [neutronNu];
    neutronDircosu = new double [neutronNu];
    neutronDircosv = new double [neutronNu];
    neutronDircosw = new double [neutronNu];
    neutronAges = new double [neutronNu];
    cmNeutronEnergies = new double [neutronNu];
    cmNeutronDircosu = new double [neutronNu];
    cmNeutronDircosv = new double [neutronNu];
    cmNeutronDircosw = new double [neutronNu];
  }
  if (photonNu > 0) { 
    photonEnergies = new double [photonNu];
    photonDircosu = new double [photonNu];
    photonDircosv = new double [photonNu];
    photonDircosw = new double [photonNu];
    photonAges = new double [photonNu];
  }
}

void cgmfEvent::allocateMemory()
{

  if (neutronNu > 0) {
    neutronEnergies = new double [neutronNu];
    neutronDircosu  = new double [neutronNu];
    neutronDircosv  = new double [neutronNu];
    neutronDircosw  = new double [neutronNu];
    neutronAges     = new double [neutronNu];

    cmNeutronEnergies = new double [neutronNu];
    cmNeutronDircosu  = new double [neutronNu];
    cmNeutronDircosv  = new double [neutronNu];
    cmNeutronDircosw  = new double [neutronNu];
  }

  if (photonNu > 0) { 
    photonEnergies = new double [photonNu];
    photonDircosu  = new double [photonNu];
    photonDircosv  = new double [photonNu];
    photonDircosw  = new double [photonNu];
    photonAges     = new double [photonNu];
  }

  if (preFissionNeutronNu>0) {
    preFissionNeutronEnergies = new double [preFissionNeutronNu];
    preFissionNeutronDircosu  = new double [preFissionNeutronNu];
    preFissionNeutronDircosv  = new double [preFissionNeutronNu];
    preFissionNeutronDircosw  = new double [preFissionNeutronNu];
    preFissionNeutronAges     = new double [preFissionNeutronNu];    
  }

}


/**********************************************************/
/*      Check Z/A/E Range                                 */
/**********************************************************/
int cgmEvent::cgmCheckRange(int z, int a, double e)
{
  if(z < 6 || z>118){
    oserr << "ERROR     :Z number " << z << " out of range\n";
    return -1;
  }
  if(a < 12 || a>300){
    oserr << "ERROR     :A number " << a << " out of range\n";
    return -1;
  }
  if(e <= 0.0){
    oserr << "ERROR     :Excitation energy " << e << " negative or zero\n";
    return -1;
  }

  double em = 0.0;
#ifdef   HAVE_PRIVATE_ENERGY_GRID
  em = custom_energy_grid[NUMBER_OF_PRIVATE_GRID-1];
#else
  em = (MAX_ENERGY_BIN-1) * ENERGY_BIN;
#endif

  if(e >= em){
    oserr << "ERROR     :Excitation energy " << e << " too high";
    return -1;
  }

  return(0);
}


cgmfEvent::cgmfEvent () {
  initialization();
  cgmAllocateMemory();
  return;
};

/************************************************************************************
             cgmfEvent: derived class constructor, destructor & methods 
 ************************************************************************************/
cgmfEvent::cgmfEvent  (int isotope, double eng, double time, double timew) : cgmEvent() {
  
	int nevents = 1;
	
  double alphaI  = 0.0;

  fissionFragmentType *lf, *hf;

  FissionFragments *ff;

  neutronNu = 0;
  photonNu  = 0;

  checkInput (isotope, eng);

  // record incident energy and ZAID of target nucleus
  incidentEnergy = eng;
  ZAIDt = isotope;

  if(firstcalc) {
    initialization();
    firstcalc = false;
    if(timew>0.0)set_time_coincidence_window(timew);
  }
 	cgmAllocateMemory();
//	allocateMemory();

  ctl.calc_montecarlo = CALC_MC;
  ctl.init_pop        = SINGLE;

  lf = new fissionFragmentType [nevents];
  hf = new fissionFragmentType [nevents];
  ff = new FissionFragments (isotope, eng, alphaI);
	ff->generateInitialFissionFragmentHistories (lf, hf, nevents);

	cout.precision(3);
	cout << std::fixed;

  double norm;
  
	//-- light fragment calc. --------------------------------------------------

  eventLF.A = lf[0].A;
  eventLF.Z = lf[0].Z;
  eventLF.U = lf[0].U;
  eventLF.spin = lf[0].spin;
  eventLF.parity = lf[0].parity;
  eventLF.KE = lf[0].KE;
  for (int i=0; i<3; i++) eventLF.preFragmentMomentum[i] = lf[0].preMomentum[i];

  specMCMainFission (&eventLF);

  //-- heavy fragment calc. --------------------------------------------------

  eventHF.A = hf[0].A;
  eventHF.Z = hf[0].Z;
  eventHF.U = hf[0].U;
  eventHF.spin   = hf[0].spin;
  eventHF.parity = hf[0].parity;
  eventHF.KE = hf[0].KE;
  for (int i=0; i<3; i++) eventHF.preFragmentMomentum[i] = hf[0].preMomentum[i];

  specMCMainFission (&eventHF);

  neutronNu = eventLF.nu + eventHF.nu;
  photonNu  = eventLF.nug + eventHF.nug;
  preFissionNeutronNu = ff->Ac0-ff->Ac;

  allocateMemory();

  if(preFissionNeutronNu>0) {  // save the momenta of the pre-fission neutrons in the lab frame
    for (int i=0; i<preFissionNeutronNu; i++) {
      double vn2=0.;
      for(int k=0;k<3;k++) vn2+=ff->labVelocities_pfn[i][k]*ff->labVelocities_pfn[i][k];
      double energy_pfn=.5*NEUTRONMASS*amuMeV*vn2;
      vn2=sqrt(vn2);
      preFissionNeutronDircosu[i] = ff->labVelocities_pfn[i][0]/vn2;
      preFissionNeutronDircosv[i] = ff->labVelocities_pfn[i][1]/vn2;
      preFissionNeutronDircosw[i] = ff->labVelocities_pfn[i][2]/vn2;
      preFissionNeutronEnergies[i] = energy_pfn;
    }
  }

	int ii;
	
	if (neutronNu>0) {
		
		for (int i=0; i<eventLF.nu; i++) {
			neutronEnergies[i]   = eventLF.neutronEnergies[i];
			neutronDircosu[i]    = eventLF.neutronDircosu[i];
			neutronDircosv[i]    = eventLF.neutronDircosv[i];
			neutronDircosw[i]    = eventLF.neutronDircosw[i];
			neutronAges[i]       = time;

      cmNeutronEnergies[i] = eventLF.cmNeutronEnergies[i];
      cmNeutronDircosu[i]  = eventLF.cmNeutronDircosu[i];
      cmNeutronDircosv[i]  = eventLF.cmNeutronDircosv[i];
      cmNeutronDircosw[i]  = eventLF.cmNeutronDircosw[i];
		}
		
		for (int i=0; i<eventHF.nu; i++) {
			ii=i+eventLF.nu;
			neutronEnergies[ii]   = eventHF.neutronEnergies[i];
			neutronDircosu[ii]    = eventHF.neutronDircosu[i];
			neutronDircosv[ii]    = eventHF.neutronDircosv[i];
			neutronDircosw[ii]    = eventHF.neutronDircosw[i];
			neutronAges[ii]       = time;

      cmNeutronEnergies[ii] = eventHF.cmNeutronEnergies[i];
      cmNeutronDircosu[ii]  = eventHF.cmNeutronDircosu[i];
      cmNeutronDircosv[ii]  = eventHF.cmNeutronDircosv[i];
      cmNeutronDircosw[ii]  = eventHF.cmNeutronDircosw[i];
		}
		
	}
	
	if (photonNu>0) {
		
		for (int i=0; i<eventLF.nug; i++) {
			photonEnergies[i]   = eventLF.photonEnergies[i];
			photonDircosu[i]    = eventLF.photonDircosu[i];
			photonDircosv[i]    = eventLF.photonDircosv[i];
			photonDircosw[i]    = eventLF.photonDircosw[i];
			photonAges[i]       = time + eventLF.photonAges[i];
		}
		
		for (int i=0; i<eventHF.nug; i++) {
			ii=i+eventLF.nug;
			photonEnergies[ii]   = eventHF.photonEnergies[i];
			photonDircosu[ii]    = eventHF.photonDircosu[i];
			photonDircosv[ii]    = eventHF.photonDircosv[i];
			photonDircosw[ii]    = eventHF.photonDircosw[i];
			photonAges[ii]       = time + eventHF.photonAges[i];
		}

	}
		
  /*** Free Allocate Memory */
  if (lf != 0) delete [] lf;
  if (hf != 0) delete [] hf;
  if (ff != 0) delete ff;


	return;
}

cgmfEvent::~cgmfEvent () {
  if (preFissionNeutronNu>0) {
    delete [] preFissionNeutronEnergies;
    delete [] preFissionNeutronDircosu;
    delete [] preFissionNeutronDircosv;
    delete [] preFissionNeutronDircosw;
    delete [] preFissionNeutronAges;
  }
  return;
}
/*****************************************************************
           end of cgmEvent constructors and destructors 
 *****************************************************************/

/*****************************************************************
 Generate fission fragment yields only.
 Initial conditions in (Z,A,KE,U,J,pi).
 *****************************************************************/
cgmfYields::cgmfYields (int isotope, double eng, int nevents, string outfilename) : cgmfEvent() {
	
	fissionFragmentType *lf, *hf;
	FissionFragments *ff;
	
	lf = new fissionFragmentType [nevents];
	hf = new fissionFragmentType [nevents];
	
	ff = new FissionFragments (isotope, eng, 0.0);
	ff->generateInitialFissionFragmentHistories (lf, hf, nevents);
	
	// save fission yields
	ofstream of;
	of.open (&outfilename[0]);
	//for (int i=0; i<nevents; i++) {
	//	of << lf[i].Z << " " << lf[i].A << " " << lf[i].KE << " " << lf[i].U << " " << lf[i].spin << " " << lf[i].parity << endl;
	//	of << hf[i].Z << " " << hf[i].A << " " << hf[i].KE << " " << hf[i].U << " " << hf[i].spin << " " << hf[i].parity << endl;
	//}

        for (int i=0; i<nevents; i++) {

         of << lf[i].Z << " " << lf[i].A << " " << lf[i].KE << " " << lf[i].U << " " << lf[i].spin << " " << lf[i].parity << " " ;

         for(int k=0;k<3;k++)

                  of << lf[i].preMomentum[k] << " " ;

         of << endl;

                of << hf[i].Z << " " << hf[i].A << " " << hf[i].KE << " " << hf[i].U << " " << hf[i].spin << " " << hf[i].parity << " " ;

                for(int k=0;k<3;k++)

                  of << hf[i].preMomentum[k] << " " ;

                of << endl;

        }
	of.close();
	
	return;
}


/*! 

*/
void cgmfEvent::checkInput (int isotope, double einc) {

  if (einc==0.0) isotope=-abs(isotope); // to ensure negative value

  int allowedIsotopes[] = {-94238,-94240,-94242,-94244,-98252,-98254,92233,92234,92235,92238,93237,94239,94241};
  int numIsotopes = sizeof(allowedIsotopes)/sizeof(int);

  bool found = false;
  for (int i=0; i<numIsotopes; i++) { if (isotope==allowedIsotopes[i]) {found=true; break;} }

  if (not found) { cerr << "ERROR: Isotope not handled by CGMF or wrong incident energy\n"; exit(-1); }
  if (isotope<0 and einc!=0.0) { cerr << "ERROR: Only spontaneous fission is handled for this isotope\n"; exit(-1); }
  if (isotope>0 and (einc>20.0 or einc<=0.0)) { cerr << "ERROR: Incident neutron energy must be greater than 0 MeV and less than or equal to 20 MeV\n"; exit (-1); }

  return;

}


/*!
 
 \brief Initialize cgmfEvent
 
 \todo Add proper documentation

 \param MAX_ENERGY_BIN: maximum number of energy bins in the continuum
 \param MAX_J: maximum number of spin values for a given energy level
 \param MAX_LEVELS: maximum number of discrete levels in a nucleus
 \param ENERG_BIN: energy bin width (in MeV)
 \param CONTINUUM_LOWER_CUT: \todo to better define!
 
 */

void cgmfEvent::initialization () {
	
	MAX_ENERGY_BIN      = 2500; //1000; // 2000;
	MAX_J               = 50; //50; // 70;
	MAX_LEVELS          = 200; //200; // 700;
	ENERGY_BIN          = 0.05; //0.1; // 0.05;
	CONTINUUM_LOWER_CUT = 0.02;
	
  outputOptions(4); // o=4
	
//  riplReadDiscreteLevelData(); // read the entire RIPL3 database of discrete levels
	readDiscreteLevelData();
  
	readkcksystdat();
	readMasses();
	readDeformations();
	readTemperatures();
	readLDP();
	readAnisotropy();
	readPreEqAngularDistributionParameters();
	readRTA();
  readSpinScaling();
	readPreeqData();

	return;
}


/*!

 \brief general memory cleanup, this needs major improvement! :MER

 */
void cgmf_cleanup(void) {
  riplDiscreteLevelsCleanup();
  cleanupRTA();
  specMCcleanup();
}

