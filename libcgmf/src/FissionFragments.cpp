/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file fissionFragments.cpp

  \brief Create Y(A,Z,TKE,U,J,pi) 

*/

/*
 * Version Log:
 *
 * v.1.2.0 (06/02/2018)
 *
 *      -- removed all XML references;
 *      -- removed option to read a main input file
 *      -- hard-coded some of the filenames and user options
 *	-- cleaned up and implemented fission yields parameterization by P. Jaffke
 */

/****** GENERAL TO DO'S */
// TODO: 1) SHOULD WE STANDARDIZE THE READ-IN FOR LEVELDENSITY, LDP, ETC. TO MATCH THINGS LIKE READLDP()?
//       2) REFORMAT THE TABLES SO ITS 1D IN ZAI AND THERE'S SOME MAP FROM ZAI --> INDEX TO CONSERVE ARRAY SIZE
//       3) GENERAL CLEANING, REMOVING UNUSED FUNCTIONS, SOME EXPLANATION IN SOME PLACES

/******************************************************************************/
/*  FISSIONFRAGMENTS.CPP                                                      */
/*        Primary function for creating Y(A,Z,TKE,J^pi) for a given fission   */
/*        reaction. Sets the yields parameterization and samples from it to   */
/*        provide a list of initial fragment conditions in masses A, charges  */
/*        Z, kinetic energies KE, excitation energies U, spins J, and         */
/*        parities PI.                                                        */
/******************************************************************************/

/* C INCLUDES */
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <functional>
// To include with XML I/O
//#include <libxml/xpath.h>
//#include <libxml/xpathInternals.h>
//#include "ffInput.hxx"

using namespace std;

/* HEADER INCLUDES */
#include "physics.h"
#include "masstable.h"
#include "config-ff.h"  // IS: added to use the MAX_NUMBER_PARAMETER_PFN parameter
#include "kinematics.h"
#include "FissionFragments.h"
//#include "Yields.h"
#include "excinterface.h"
#include "evapinterface.h"
#include "terminate.h"

// from CGM
#include "config.h"
#include "cgm.h"
#include "ripl2levels.h"
//#include "global_var.h"

extern std::function< double(void) > rng_cgm;

// Default constructor
FissionFragments::FissionFragments(void) { }

/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Given a nucleus ZAID = Z*1000 + A and an incident neutron energy Einc (MeV) 
 * and an (optional) user-specified spin-scaling factor alphaSpin initialize
 * set up the fission fragments for this fission reaction.
 ******************************************************************************/
FissionFragments::FissionFragments(int ZAID, double Einc, double alphaSpin) {

    barrier = NULL;
    sp_pfn_mult = NULL;
    temperatures = NULL;
    levelDensities = NULL;
    levelDensityParameters = NULL;
   
    incidentEnergy  = Einc; // Incident energy (which may be changed by multi-chance fission)
    incidentEnergy0 = Einc; // Initial incident energy Einc (doesn't change)

    // Target nucleus
    ZAIDt = ZAID;
    Zt = int(ZAIDt/1000.0);
    At = ZAIDt - 1000*Zt;
    alphaI0 = 0.0;

    // Compound (parent/fissioning) nucleus
    if (incidentEnergy==0.0) { // spontaneous fission
	Zc = Zt; Ac = At;
    } else { // neutron-induced fission
	Zc = Zt; Ac = At+1;
    }
    ZAIDc  = 1000*Zc + Ac;
    ZAIDc0 = ZAIDc;	// Initial compound nucleus (doesn't change)
    Zc0 = Zc;		// Initial compound charge (doesn't change)
    Ac0 = Ac;		// Initial compound mass (doesn't change)
   
    // Defines the symmetric point in A/Z
    Asym = Ac/2; Zsym = Zc/2;
   
    // Masses and neutron separation energies
    massExcessCompound = getMassExcess(ZAIDc);
    SnCompound = 0.0; // Spontaneous fission
    if (incidentEnergy>0.0) SnCompound = -getMassExcess(ZAIDc) + getMassExcess(ZAIDc-1) + getMassExcess(1);

    // Sets the options given the ZAID
    setOptions();

    // Initializes the yields
    init();

    // If a user-specified alphaSpin was provided, use this for alphaI0
    if (alphaSpin>0) alphaI0 = alphaSpin;
    return;
}

/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Constructs an instance of FissionFragments based on user input:
 * ZAID:       ZAID number of target nucleus
 * Einc:       incident neutron energy in MeV
 * alphaSpin:  spin distribution factor
 * yieldsFile: file providing the initial fission fragment yields in (A,Z,TKE,Y)
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAID, double Einc, double *alphaSpin, string yieldsFile) {

    incidentEnergy = Einc;

    // Target nucleus
    ZAIDt = ZAID;
    Zt = int(ZAIDt/1000.0);
    At = ZAIDt-1000*Zt;

    //-- compound (parent, fissioning) nucleus
    if (incidentEnergy==0.0) { // spontaneous fission
	Zc = Zt; Ac = At;
    } else { // neutron-induced fission
	Zc = Zt; Ac = At + 1;
    }
    ZAIDc = 1000*Zc + Ac;

    //-- symmetric fragment
    Asym = Ac/2; Zsym = Zc/2;

    // setting options

    yieldsOption = "YAZTKEfile";
    dZ = 1;
    sigmaZ = 0;
    alphaI =* alphaSpin;
  
    energySortingOption = "RTA";
  
    yieldsFilename = yieldsFile;
  
    init();
    return;
}




/*******************************************************************************
 * Constructor
 *------------------------------------------------------------------------------
 * Define an input filename that corresponds to a ZAID and incidentEnergy. Then
 * reads the corresponding XML input file and builds the FF yields.
 ******************************************************************************/
FissionFragments::FissionFragments (int ZAIDf, double excitationEnergy, string option)
{
    if (option.compare("single")!=0) { cerr << "Wrong option for FissionFragments\n"; }
    // no need to do anything
    return;
}

/*******************************************************************************
 * setOptions
 *------------------------------------------------------------------------------
 * Sets the options for the fission fragment parameterizations using the target
 * ZAIDt and incidentEnergy. Determines what type of systematics (yieldsOption)
 * should be used, a yieldsFilename (if yieldsOption="YATKE"), the number of
 * charges dZ to be considered around the most-probable charge and the
 * energySortingOption ("RTA" for R_T(A) and "fixedRT" for a constant value).
 * We note that dZ and dZ1 are used to account for the change in Zp and sigmaZ
 * with increasing incident neutron energy. Thus, dZ1 is used to define the
 * level-density parameter, temp., etc. tables and dZ to define the Y(Z|A).
 ******************************************************************************/
void FissionFragments::setOptions(void) {
   
    yieldsFilename = datadir; // To prepend to the yieldsFilename

    sf_flag = false;
    anisotropy_flag = false;
    recomputeYields = true; // Default is to construct the yields (instead of reading them)

    switch (ZAIDt) {

	/* URANIUMS */
	case 92233: // U233
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    break;

	case 92234: // U234
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "fixedRT";
	    anisotropy_flag = true;
	    RT = 1.2;
	    break;

	case 92235: // U235
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    break;

	case 92238: // U238
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ+3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    if (incidentEnergy==0.0) {
		cerr<<"ERROR: Only neutron-induced fission allowed for U238!"<<endl; exit(0);
	    }
	    break;

	/* NEPTUNIUMS */
	case 93237: // Np237
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    break;

	/* PLUTONIUMS */
	case 94238: // Pu238
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Pu238!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 6;
	    dZ1 = dZ;
	    energySortingOption = "fixedRT";
	    RT = 1.2;
	    break;

	case 94239: // Pu239
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    break;

	case 94240: // Pu240
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Pu240!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 6;
	    dZ1 = dZ;
	    energySortingOption = "RTA";
	    break;

	case 94241: // Pu241
	    yieldsOption = "Systematics";
	    dZ  = 3;
	    dZ1 = dZ + 3;
	    energySortingOption = "RTA";
	    anisotropy_flag = true;
	    break;

	case 94242: // Pu242
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Pu242!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 6;
	    dZ1 = dZ;
	    energySortingOption = "RTA";
	    break;

	case 94244: // Pu244
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Pu244!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 6;
	    dZ1 = dZ;
	    energySortingOption = "fixedRT";
	    RT = 1.2;
	    break;

	/* CALIFORNIUMS */
	case 98252: // Cf252
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Cf252!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 7;
	    dZ1 = dZ;
	    energySortingOption = "RTA";
	    break;

	case 98254: // Cf254
	    if (incidentEnergy!=0.0) {
		cerr<<"ERROR: Only spontaneous fission allowed for Cf254!"<<endl; exit(0);
	    } else {
		sf_flag = true;
	    }
	    yieldsOption = "Systematics";
	    dZ  = 6;
	    dZ1 = dZ;
	    energySortingOption = "RTA";
	    break;

	/* NOT SUPPORTED */
	default:
	    cerr << "ERROR: Cannot use CGMF (yet) for this isotope and/or incident energy" << endl;
	    exit(0);
    }
    return;
}

/*******************************************************************************
 * Init
 *------------------------------------------------------------------------------
 * First main routine to be called after constructor
 ******************************************************************************/
void FissionFragments::init () {

  resetYields ();
  buildYields ();
  readAnisotropy ();
	
}

/*
	 readTemperatures
*/
// void FissionFragments::readTemperatures () {
	
// 	ifstream fp;

// 	string str = TEMPERATURESFILE;

// 	cout << str << endl;

// 	fp.open(&str[0]); // try current directory
// 	if (!fp) { // then system data directory
// 		str = datadir + str;
// 		fp.open(&str[0]);
// 	}

// 	if (!fp) cgmTerminateCode("Temperature data file not found");

//     while(getline(fp,str)) {
// 		if(str[0] == '#') continue;
// 		int c=-1;
// 		while (!fp.eof()) {
// 			c++;
// 			fp >> ZAIDlist[c];
// 			for (int i=0; i<32; i++) fp >> temp[c][i];
// 		}
// 	}

// 	return;
	
// }

	 
/*******************************************************************************
 * Set yields YA[], YZ[], ... to zero.
 ******************************************************************************/
void FissionFragments::resetYields() {
  
  for (int i=0; i<NUMA; i++) {
    YA[i]=0.0;
    maxYieldZ[i]=0.0;
    QfA[i]=0.0;
  }
  
  for (int i=0; i<NUMZ; i++) YZ[i]=0.0;
  
  for (int i=0; i<NUMTKE; i++) {
    YTKE[i]=0.0;
    QfTKE[i]=0.0;
  }
  
  for (int i=0; i<NUMA; i++) {
    for (int j=0; j<NUMTKE; j++) {
      YATKE[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMdZ; i++) {
    for (int j=0; j<NUMA; j++) {
      YZA[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMZ; i++) {
    for (int j=0; j<NUMA; j++) {
      YZA2[i][j] = 0.0;
    }
  }
  
  for (int i=0; i<NUMZ; i++) {
    for (int j=0; j<NUMA; j++) {
      for (int k=0; k<NUMTKE; k++) {
        YZATKE[i][j][k] = 0.0;
      }
    }
  }  
  
}

/*******************************************************************************
 * Reads a main input file without XML formatting.
 ******************************************************************************/
void FissionFragments::readMainInputFile (string inputFilename) {
  
  ifstream inputFile;
  string line, key, value;
  size_t found;
  
  inputFile.open(&inputFilename[0]);
  if (!inputFile) { cerr << "[readMainInputFile] Cannot find main input file!\n"; cerr << inputFilename << "\n"; exit(-1); }
  
  while (!inputFile.eof()) {
    
    inputFile >> line;
    
    found = line.find("=");
    if (found!=string::npos) {
      
      key = line.substr(0, int(found));
      value = line.substr(int(found)+1,80);
      //cout << key << " : " << value << "\n";
      
      if (key == "At") At = atoi(value.c_str());
      if (key == "Zt") Zt = atoi(value.c_str());
      if (key == "Ap") Ap = atoi(value.c_str());
      if (key == "Zp") Zp = atoi(value.c_str());
      
      if (key == "incidentEnergy") incidentEnergy = atof(value.c_str());
      
      if (key == "yieldsOption") yieldsOption = value.c_str();
      if (key == "yieldsFilename") yieldsFilename = value.c_str();
      
      //      if (key == "TKEmin") TKEmin = atoi(value.c_str());
      //      if (key == "TKEmax") TKEmax = atoi(value.c_str());
      
      if (key == "Amin") Amin = atoi(value.c_str());
      if (key == "Amax") Amax = atoi(value.c_str());
      if (key == "Zmin") Zmin = atoi(value.c_str());
      if (key == "Zmax") Zmax = atoi(value.c_str());
      
      if (key == "dZ") dZ = atoi(value.c_str());
      if (key == "sigmaZ") sigmaZ = atof(value.c_str());
      
      if (key == "energySortingOption") energySortingOption = value;
      
      if (key == "RT") RT = atof(value.c_str());
      
      if (key == "alphaI") alphaI=atof(value.c_str());
      
    }
    
  }
  
  
  
  //-- compound (parent, fissioning) nucleus
  if (incidentEnergy==0.0) { // spontaneous fission
    Zc=Zt; Ac=At;
  } else { // neutron-induced fission
    Zc=Zt; Ac=At+1;
  }
  
  //-- symmetric fragment
  Asym=Ac/2; Zsym=Zc/2;
  
  inputFile.close();
  
}

// TO INCLUDE WITH XML I/O
/*******************************************************************************
 * First reads an XML input data file, and creates an instance of FissionFragments.
 * It then builds the fission fragment yields Y(A,Z,TKE).
 *
 * << XCode note >>
 *
 * It requires XERCES and XSD C++ libraries installed and loaded as "external
 * frameworks and libraries". Also, paths to the libraries are to be included
 * in the project header files path search.
 *
 ******************************************************************************/
/* void FissionFragments::readXMLInputFile (string inputFilename) {
 
 try {
 
 auto_ptr<ffInputType> mainInput (ffInput(inputFilename)); // parses XML data file using XSD C++ routines
 
 //-- target nucleus
 Zt = mainInput->target().Z();
 At = mainInput->target().A();
 
 incidentEnergy = mainInput->incidentEnergy();
 
 //-- compound (parent, fissioning) nucleus
 if (incidentEnergy==0.0) { // spontaneous fission
 Zc=Zt; Ac=At;
 } else { // neutron-induced fission
 Zc=Zt; Ac=At+1;
 }
 
 //-- symmetric fragment
 Asym=Ac/2; Zsym=Zc/2;
 
 TKEmin = mainInput->yields().TKEmin();
 TKEmax = mainInput->yields().TKEmax();
 
 Amin = mainInput->yields().Amin();
 Amax = mainInput->yields().Amax();
 
 Zmin = mainInput->yields().Zmin();
 Zmax = mainInput->yields().Zmax();
 
 dZ     = mainInput->yields().dZ();
 sigmaZ = mainInput->yields().sigmaZ();
 // dZFile = mainInput->yields().dZFile();
 
 energySortingOption = mainInput->energySorting().option();
 RT = mainInput->energySorting().RT();
 
 yieldsOption   = mainInput->yields().yieldsOption();
 yieldsFilename = mainInput->yields().yieldsFilename();
 
 }
 
 catch (const xml_schema::exception& e)
 {
 cerr << e << endl;
 exit(1);
 }
 
 
 }
 */

// Destructor
FissionFragments::~FissionFragments () {

//  if (fragments != 0) delete [] fragments;
  if (n_pfn > 0) {
    if (sp1 != 0) delete [] sp1;
    if (sp != 0) delete [] sp;
  }
  if (tc != 0) {
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[k].tran;
    delete [] tc;
  }
  if (barrier != 0) delete [] barrier;
  if (emissprob != 0) delete [] emissprob;
  if (energy_grid1 != 0) delete [] energy_grid1;
  if (sp_pfn_mult != 0) {
    for (int i=0; i<nemit; i++) {
      if (sp_pfn_mult[i] != 0) delete [] sp_pfn_mult[i];
    }
    delete [] sp_pfn_mult;
  }
  if (temperatures != 0) {
    for (int i=Amin; i<=Amax; i++) {
      if (temperatures[i] != 0) {
        for (int j=0; j<2*dZ; j++) {
          if (temperatures[i][j] != 0) delete [] temperatures[i][j];
        }
        delete [] temperatures[i];
      }
    }
    delete [] temperatures;
  }
  if (levelDensities != 0) {
    for (int i=Amin; i<=Amax; i++) {
      if (levelDensities[i] != 0) {
        for (int j=-dZ; j<=dZ; j++) {
          if (levelDensities[i][j+dZ] != 0) delete [] levelDensities[i][j+dZ];
        }
        delete [] levelDensities[i];
      }
    }
    delete [] levelDensities;
  }
  if (levelDensityParameters != 0) {
    for (int i=Amin; i<=Amax; i++) {
      if (levelDensityParameters[i] != 0) {
        for (int j=-dZ; j<=dZ; j++) {
          if (levelDensityParameters[i][j+dZ] != 0) delete [] levelDensityParameters[i][j+dZ];
        }
        delete [] levelDensityParameters[i];
      }
    }
    delete [] levelDensityParameters;
  }

}

/*!
 Reads deformation energies at scission calculated by P.Moller, 2012.
 */

void FissionFragments::readMollerDeformationEnergies (string filename) {
  
  int A;
  double Edef;
  
  std::fill_n(EdefA, NUMA, 0.0);
  
  ifstream data;
  data.open (&filename[0]);
  
  while (!data.eof()) {
    data >> A >> Edef;
    EdefA[A] = Edef;
  }
  data.close();
  
}

// TODO: REFORMAT LEVELDENSITIES TO BE A 1D ARRAY OF ZAID? WITH SOME MAPPING TO CONSERVE SPACE? SIMILAR TO readTemperatures. (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * readLevelDensityTables
 ******************************************************************************/
void FissionFragments::readLevelDensityTables (string datafile) {
  
  cout << "[readLevelDensityTables] Running... ";
  
  
  ifstream data;
  string str;
  unsigned pos;
	int A;
//  int Z, deltaZ;
  
  levelDensities = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    levelDensities[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      levelDensities[i][j+dZ] = new double [NUME];
      std::fill_n(levelDensities[i][j+dZ], NUME, 0.0);
    }
  }
  
  // read data file
  data.open (&datafile[0]);
  while (getline(data,str)) {
    
    if (str.find("#")==0) {
      
      pos = (int) str.find("A=");  A  = atoi(str.substr(pos+2,3).c_str());
//      pos = str.find("Zp="); Z  = atoi(str.substr(pos+3,3).c_str());
//      pos = str.find("dZ="); deltaZ = atoi(str.substr(pos+3,2).c_str());
      
      if (A>Amax || A<Amin) { continue; }
      
      //			cout << A << " " << Z << " " << deltaZ << "\n";
      
      getline(data,str); getline(data,str); // skip two lines
      
      getline(data,str);
      while (!str.empty()) {
        for (int k=0; k<NUME; k++) {
          //					data >> levelDensityParameters[A][0][k];
          //					cout << k << " " << levelDensityParameters[A][0][k] << "\n";
          for (int j=-dZ; j<=dZ; j++) {
            levelDensities[A][j+dZ][k] = atof(str.substr((j+dZ+1)*10,10).c_str());
            //						cout << setprecision(6) << atof(str.substr((j+dZ+1)*10,11).c_str()) << "\n";
          }
          getline(data,str);
        }
      }
      
      // TO BE CONTINUED... NEED TO SOLVE PROBLEMS OF DIFFERENT ENERGY GRID
      // AND DIFFERENT FRAGMENTS PRODUCED IN DIFFERENT REACTIONS
      
    }
  }
  data.close();
  
  cout << "OK\n";
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: REFORMAT LEVELDENSITIES TO BE A 1D ARRAY OF ZAID? WITH SOME MAPPING TO CONSERVE SPACE? SIMILAR TO readTemperatures. (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * computeLevelDensityTables
 *------------------------------------------------------------------------------
 * Computes the level density tables rho(A,Z,U) for all of the fragment (A,Z)
 * pairs and over a given energy grid (with NUME energy points). The table
 * rho(A,Z,U) is used in the excitation energy sorting.
 ******************************************************************************/
void FissionFragments::computeLevelDensityTables() {
   
    Nucleus nucleus;
    ofstream out;
    string filename;
    stringstream line;

    /* FOR SAVING THE LEVEL DENSITY TABLE COMMENT OUT THE LINES BELOW

    out.open("ldtables.dat", ios::out);
    out << setprecision(3) << setiosflags(ios::scientific);
   */

    cout << "[computeLevelDensityTables] Running... ";

    levelDensities = new double **[NUMA];
    for (int i=Amin; i<=Amax; i++) {
	levelDensities[i] = new double *[NUMdZ];
	for (int j=-dZ1; j<=dZ1; j++) {
	    levelDensities[i][j+dZ1] = new double [NUME];
	    nucleus.za.setZA((int) floor(Z0[i]+j+0.5),i);
	    nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign); // Use reassign
	    statFixDiscreteLevels(&nucleus);
	    statSetupEnergyBin(&nucleus);
	    statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
	    for (int k=0; k<NUME; k++) {
		levelDensities[i][j+dZ1][k] = ldLevelDensity(ldEnergyGrid[k], i, &nucleus.ldp);
	    }
	}

      /* FOR SAVING THE TABLE
	out << "\n# A=" << i << " ; Zp=" << int(Z0[i]+0.5) << " ; dZ=" << dZ << "\n";
	out << "# U / Z = ";
	out << setw(10);
	for (int j=-dZ; j<=dZ; j++) { out << setw(10) << int(Z0[i]+0.5)+j; }
	out << endl << endl;
	for (int k=0; k<NUME; k++) {
	    out << setw(10) << ldEnergyGrid[k];
	    for (int j=-dZ; j<=dZ; j++) {
		out << setw(10) << levelDensities[i][j+dZ][k];
	    }
	out << "\n";
	}
      */
    }
    //   out.close();
    cout << "OK\n";
    return;   
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: DO WE NEED THIS ANYMORE SINCE WE HAVE READLDP()? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * readLevelDensityParameterTables
 ******************************************************************************/
void FissionFragments::readLevelDensityParameterTables (string datafile) {
  
  cout << "[readLevelDensityParameterTables] Running... ";
  
  
  ifstream data;
  string str;
  unsigned pos;
	int A;
//  int Z, deltaZ;
  
  levelDensityParameters = new double **[NUMA];
  for (int i=Amin; i<=Amax; i++) {
    levelDensityParameters[i] = new double *[NUMdZ];
    for (int j=-dZ; j<=dZ; j++) {
      levelDensityParameters[i][j+dZ] = new double [NUME];
      std::fill_n(levelDensityParameters[i][j+dZ], NUME, 0.0);
    }
  }
  
  // read data file
  data.open (&datafile[0]);
  while (getline(data,str)) {
    
    if (str.find("#")==0) {
      
      pos = (int) str.find("A=");  A  = atoi(str.substr(pos+2,3).c_str());
//      pos = (int) str.find("Zp="); Z  = atoi(str.substr(pos+3,3).c_str());
//      pos = (int) str.find("dZ="); deltaZ = atoi(str.substr(pos+3,2).c_str());
      
      if (A>Amax || A<Amin) { continue; }
      
      getline(data,str); getline(data,str); // skip two lines
      
      getline(data,str);
      while (!str.empty()) {
        for (int k=0; k<NUME; k++) {
          //					data >> levelDensityParameters[A][0][k];
          //					cout << k << " " << levelDensityParameters[A][0][k] << "\n";
          for (int j=-dZ; j<=dZ; j++) {
            levelDensityParameters[A][j+dZ][k] = atof(str.substr((j+dZ+1)*10,10).c_str());
            //						cout << setprecision(6) << atof(str.substr((j+dZ+1)*10,11).c_str()) << "\n";
          }
          getline(data,str);
          //					cout << str << "\n";
        }
      }
      
      // TO BE CONTINUED... NEED TO SOLVE PROBLEMS OF DIFFERENT ENERGY GRID
      // AND DIFFERENT FRAGMENTS PRODUCED IN DIFFERENT REACTIONS
      
    }
  }
  data.close();
  
  cout << "OK\n";
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: DO WE NEED THIS ANYMORE SINCE WE HAVE READLDP()? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * computeLevelDensityParameterTables
 *------------------------------------------------------------------------------
 * Computes level density parameter tables f(A,Z,U), which are later used to
 * compute energy sorting at scission, if energySortingOption = "fixedRT".
 ******************************************************************************/
void FissionFragments::computeLevelDensityParameterTables() {
   
    Nucleus nucleus;
    ofstream out;
    string filename;
    stringstream line;

    /* FOR SAVING THE LEVEL DENSITY PARAMETER TABLE

    out.open("ldparamtables.dat", ios::out);
    out << setprecision(6) << setiosflags(ios::fixed);
    */

    cout << "[computeLevelDensityParameterTables] Running... ";

    levelDensityParameters = new double **[NUMA];
    for (int i=Amin; i<=Amax; i++) {

	levelDensityParameters[i] = new double *[NUMdZ];
	for (int j=-dZ1; j<=dZ1; j++) {
	    levelDensityParameters[i][j+dZ1] = new double [NUME];
	    nucleus.za.setZA((int)floor(Z0[i])+j,i);
	    nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign); // Use reassign
	    statFixDiscreteLevels(&nucleus);
	    statSetupEnergyBin(&nucleus);
	    statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);
	    for (int k=0; k<NUME; k++) {
		levelDensityParameters[i][j+dZ1][k] = ldDensityParameter(ldEnergyGrid[k], i, &nucleus.ldp);
	    }
	}

      /* FOR SAVING THE TABLE
	out << "\n# A=" << i << " ; Zp=" << int(Z0[i]+0.5) << " ; dZ=" << dZ << "\n";
	out << "# U / Z = ";
	out << setw(10);
	for (int j=-dZ; j<=dZ; j++) { out << setw(10) << int(Z0[i]+0.5)+j; }
	out << endl << endl;
	for (int k=0; k<NUME; k++) {
	    out << setw(10) << ldEnergyGrid[k];
	    for (int j=-dZ; j<=dZ; j++) {
		out << setw(10) << levelDensityParameters[i][j+dZ][k];
	    }
	    out << "\n";
	}
      */
    }
    //   out.close();
    cout << "OK" << endl;
    return;
}
// -----------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************
 * computeEnergyFromMaxEntropy
 *------------------------------------------------------------------------------
 * Compute mean light fragment <Ul> for a given total intrinsic excitation
 * energy, following KHS formula [Phys. Rev. C83, 061601(R) (2011)].
 ******************************************************************************/
double FissionFragments::computeEnergyFromMaxEntropy(int Zl, int Al, double TIXE) {
    // TIXE: Total Intrinsic eXcitation Energy (MeV)
    // (Zl,Al): mass and charge of light fragment

    double Ul;
    int Ah, Zh;

    Ah = Ac - Al;
    Zh = Zc - Zl;

    int dZl, dZh;

    dZl = int ( Zl - floor(Z0[Al] + 0.5) + dZ1 );
    dZh = int ( Zh - floor(Z0[Ah] + 0.5) + dZ1 );

    int kl, kh;
    double numerator, denominator;

    int maxIndex = (int) (TIXE/deltaE) + 1;

    numerator = 0.0;
    denominator = 0.0;
    for (kl=1; kl<=maxIndex; kl++) {
	Ul = (kl-1)*deltaE; // deltaE=0.25 (see config-ff.h)
	kh =  (int) (floor((TIXE-Ul)/deltaE) + 1);
	numerator += Ul * levelDensities[Al][dZl][kl] * levelDensities[Ah][dZh][kh]; // TODO: Change to 1D array in ZAI? (P.J.)
	denominator += levelDensities[Al][dZl][kl] * levelDensities[Ah][dZh][kh];
    }
    if (denominator!=0) {
	return numerator/denominator;
    } else {
	cout << "[computeEnergyFromMaxEntropy] denominator is null!\n";
	cout << Al << " " << Zl << "\n";
	return 0;
    }
}

/*******************************************************************************
 * studyEnergySorting()
 *------------------------------------------------------------------------------
 ******************************************************************************/
void FissionFragments::studyEnergySorting () {
  int i, count;
  string ldFilename, tempFilename, ctFilename;
  
  ldFilename   = datadir + "levelDensities.CGMF.out";
  tempFilename = datadir + "temperatures.CGMF.out";
  ctFilename   = datadir + "constantTemperatures.CGMF.out";
		
  // LEVEL DENSITIES
  ofstream ldFile;
  ldFile.open (&ldFilename[0]);
  count=-1;
  for (i=Amin; i<Amax; i++) {
    for (int j=-dZ; j<=dZ; j++) {
      ldFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << (int) floor(Z0[i]+j+0.5) << " (j=" << j << ")\n\n";
      for (int k=0; k<NUME; k++) {
        ldFile << ldEnergyGrid[k] << " " << levelDensities[i][j+dZ][k] << " " << log(levelDensities[i][j+dZ][k]) << "\n"; // TODO: Change to 1D array in ZAI? (P.J.)
      }
      ldFile << "\n";
    }
  }
  ldFile.close();
  
  
  // TEMPERATURES
  ofstream tempFile;
  tempFile.open (&tempFilename[0]);
  count=-1;
  for (int i=Amin; i<=Amax; i++) {
    for (int j=0; j<=2*dZ; j++) {
      tempFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << (int) floor(Z0[i]+j-dZ+0.5) << "\n\n";
      tempFile << setprecision(3) << setiosflags(ios::scientific);
      for (int k=0; k<NUME; k++) {
        tempFile << ldEnergyGrid[k] << " " << temperatures[i][j][k] << " " << levelDensities[i][j][k] << "\n";
      }
      tempFile << "\n";
    }
  }
  tempFile.close();
  
  // LOW-ENERGY TEMPERATURES AS A FUNCTION OF MASS
  ofstream ctFile;
  ctFile.open (&ctFilename[0])  ;
  double t, sum;
  for (int i=Amin; i<=Amax; i++) {
    t=0.0;
    sum=0.0;
    for (int j=0; j<=2*dZ; j++) {
      t += temperatures[i][j][0] * YZA[j][i];
      sum += YZA[j][i];
    }
    if (sum!=0.0) t /= sum;
    ctFile << i << " " << t << "\n";
  }
  ctFile.close();
  
  ofstream maxEntropyFile;
  maxEntropyFile.open ("maxEntropy.CGMF.out");
  // Compute <Ul> from max. entropy principle
  int A1;
  int Z1;
  
  double x, TXE;
  count=-1;
  for (int i=70; i<=Asym; i++) {
    A1 = i;
    Z1 = (int) floor(Z0[i]+0.5);
    maxEntropyFile << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << Z1 << "\n\n";
    //	  cout << "# [gnuplot #" << ++count << "] A=" << i << " ; Z=" << Z1 << "\n\n";
    for (int i=1; i<=200; i++) {
      TXE = (i-1)*deltaE;
      x = computeEnergyFromMaxEntropy(Z1, A1, TXE);
      maxEntropyFile << TXE << " " << x << " " << TXE/(1+float(Ac-A1)/float(A1)) << "\n";
    }
    maxEntropyFile << "\n";
  }
  
  maxEntropyFile.close();
  
  cout << "TMP\n";
  exit(0);
  
}

/*******************************************************************************
 * generateInitialFissionFragmentHistories
 *------------------------------------------------------------------------------
 * Generates the Monte Carlo histories of primary fission fragments in mass A,
 * charge Z, kinetic energy KE, excitation energy U, spin J, and parity PI.
 * Does this for numberEvents events and stores in lightFragment, heavyFragment
 ******************************************************************************/
void FissionFragments::generateInitialFissionFragmentHistories(fissionFragmentType* lightFragment, fissionFragmentType* heavyFragment, const int numberEvents) {

  // calculate the anisotropy and CDF
  computeFragmentAngularDistribution();

  int i = 0; // Keeps track of the event number

  if (yieldsOption == "Systematics") {
   	Yields ffy(ZAIDc,sf_flag);
    while (i<numberEvents) {
      index_save = i;
      if (recomputeYields) constructYields(&ffy);
      sampleFissionFragments(lightFragment+i,heavyFragment+i,&ffy);
      if (lightFragment[i].U != 0.0 && heavyFragment[i].U != 0.0) ++i;
    }
  } else {
    while (i<numberEvents) {
	    index_save = i;
	    sampleFissionFragments(lightFragment+i,heavyFragment+i);
	    if (!(lightFragment[i].U == 0.0 && heavyFragment[i].U == 0.0)) i++;
    }
  }
  return;
}

/*******************************************************************************
 * generateSingleFissionFragments
 *------------------------------------------------------------------------------
 * Used for studying the decay for a specific fragment (ZAIDf) at an excitation
 * energy Eexc. Stores it and its complement in lightFragment, heavyFragment
 * for numberEvents events.
 ******************************************************************************/
void FissionFragments::generateSingleFissionFragments(fissionFragmentType* lightFragment, fissionFragmentType* heavyFragment, int ZAIDf, double Eexc, const int numberEvents) {

    int Zf = int(ZAIDf/1000);
    int Af = ZAIDf%1000;

    for (int i=0; i<numberEvents; i++) {
	lightFragment[i].A = Af;
	lightFragment[i].Z = Zf;
	lightFragment[i].U = Eexc;
	lightFragment[i].KE = 0.0;
	lightFragment[i].spin = 0;
	lightFragment[i].parity = 0;
	heavyFragment[i].A = Af;
	heavyFragment[i].Z = Zf;
	heavyFragment[i].U = Eexc;
	heavyFragment[i].KE = 0.0;
	heavyFragment[i].spin = 0;
	heavyFragment[i].parity = 0;
    }
    return;
}


/*******************************************************************************
 * generateInitialFissionFragmentHistories (file out)
 *------------------------------------------------------------------------------
 * Generates Monte Carlo histories of primary fission fragment with their
 * initial excitation energy and intial spin. The number of generated histories
 * is a user input (numberEvents).
 ******************************************************************************/
void FissionFragments::generateInitialFissionFragmentHistories (string outputFilename, int numberEvents) {
  fissionFragmentType lightFragment, heavyFragment;
	Yields *ffy;
	
	ofstream outputFile;
  outputFilename = datadir+outputFilename;
  outputFile.open(&outputFilename[0]);
  outputFile.precision(2);
  outputFile << fixed;
  computeFragmentAngularDistribution(); // compute the anisotropy distribution
	
	if (yieldsOption == "Systematics") { // Use of fission fragment yield systematics

		if (incidentEnergy==0.0) {
			ffy = new Yields(ZAIDc, true);
		} else {
			ffy = new Yields(ZAIDc, false);
		}
		for (int i=1; i<=numberEvents; i++) { // << 1.0.6 >>
		
			index_save=i-1;
			constructYields (ffy);
			sampleFissionFragments (&lightFragment, &heavyFragment, ffy);
			
			if (lightFragment.U>0.0 && heavyFragment.U>0.0) {
			
				if (i>1) outputFile << endl;
			
				outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
				setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
				lightFragment.spin << setw(3) << lightFragment.parity << endl;
			
				outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
				setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
				heavyFragment.spin << setw(3) << heavyFragment.parity;
			}
		
		}
		
	} else { // use of specific fission fragment yields
		
		for (int i=1; i<=numberEvents; i++) {
			
      sampleFissionFragments (&lightFragment, &heavyFragment);
      
			if (lightFragment.U>0.0 && heavyFragment.U>0.0) {
      
				if (i>1) outputFile << endl;
      
				outputFile << setw(4) << lightFragment.Z << setw(4) << lightFragment.A <<
				setw(10) << lightFragment.KE << setw(8) << lightFragment.U << setw(8) <<
				lightFragment.spin << setw(3) << lightFragment.parity << endl;
      
				outputFile << setw(4) << heavyFragment.Z << setw(4) << heavyFragment.A <<
				setw(10) << heavyFragment.KE << setw(8) << heavyFragment.U << setw(8) <<
				heavyFragment.spin << setw(3) << heavyFragment.parity;
			}
    
		}
  
	}
	
  outputFile.close();
  
}

/*******************************************************************************
 * sampleFissionFragments
 *------------------------------------------------------------------------------
 * Performs the Monte Carlo sampling of the fission fragment yields Y(A,Z,TKE)
 * to select a fragment pair (Al,Zl,KEl) - (Ah,Zh,KEh) and then shares the
 * excitation energy and select the J^Pi for each fragment via the function
 * "computeFragmentInitialConditions". Assumes Y(A,TKE) and Y(Z|A) have been
 * built already.
 ******************************************************************************/
void FissionFragments::sampleFissionFragments(fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment) {

    // Fragment properties
    double TKE;
    double maxYz;
    double KEl, KEh;
    int iTKE, iA, iZ; iTKE = iA = iZ = 0; // Index version of TKE, A, Z
    int j;
    int Al, Zl;
    int Ah, Zh;
    int iZmin, iZmax; iZmin = -dZ; iZmax = dZ;

    // Excitation energies, spins, and parities
    double Ul, Uh;
    double spinl, spinh;
    int parl, parh;

    double y, yz; // Yields

    double energyRelease; // Q-value for a given fragment pair

    double *momentum; momentum = new double [3]; // (x,y,z) components of fragment momentum

    int i = 0; // Counter for the number of (A,TKE) sampling attempts    

    // Monte Carlo sampling
    do {
	i++;

	// First select the TKE and A
	TKE  = rng_cgm()*(TKEmax - TKEmin + 1) + TKEmin;
	iTKE = int(floor(TKE));
	iA   = int(floor(rng_cgm()*(Amax - Amin + 1) + Amin));

	y = rng_cgm()*maxYieldATKE; // Sample a yield

	// Check if the sampled yield is lower than the Y(A,TKE) for the sampled A and TKE
	if (y <= YATKE[iA][iTKE]) {

	    // If our sampling was successful, then proceed to get the Z
	    maxYz = 1.01*maxYieldZ[iA];
	    j = 0; // Counter for the number of Z sampling attempts
	    do {
		j++;
		// Sample a Z value within dZ
		iZ = int(floor(rng_cgm()*(iZmax - iZmin + 1) + iZmin));
		yz = rng_cgm()*maxYz; // Sample a yield

		// Don't sample over 10000 times for a single Z
		if (j>=10000) { cerr << "Charge sampling exceeded!" << endl; }
	    } while (yz > YZA[iZ + dZ][iA]);
	}

	// Don't sample over 10000 times for a single event
	if (i>=10000) { cerr << "[sampleFissionFragments] Sampling did not converge!" << endl; }
    } while (y > YATKE[iA][iTKE]);

    // TODO: CAREFUL WE CAN SAMPLE AN A<ASYM AND A Z>ZSYM, SO WHICH FRAGMENT IS THIS? (P.J.)
    // -----------------------------------------------------------------------------------------------------------------------------
    // Check if we've sampled a light or heavy fragment
    if (iA<=Ac/2) {
	Al = iA;
	Zl = iZ + int(floor(Z0[iA] + 0.5));
    } else {
	Al = Ac - iA;
	Zl = Zc - iZ - int(floor(Z0[iA] + 0.5));
    }
    Zh = Zc - Zl;
    Ah = Ac - Al;
    // -----------------------------------------------------------------------------------------------------------------------------

    // Compute the excitation energy sharing and the spin/parity assignment
    computeFragmentInitialConditions(TKE,Zl,Al,Zh,Ah,&Ul,&Uh,&spinl,&spinh,&parl,&parh,&energyRelease,lightFragment->preMomentum,heavyFragment->preMomentum);

    for (unsigned int i=0;i<3;i++) { // Initialize for now
	lightFragment->postMomentum[i] = lightFragment->preMomentum[i];
	heavyFragment->postMomentum[i] = heavyFragment->preMomentum[i];
    }

    // ground-state masses of the fragments << PT >>
    double Ml = Al*amuMeV + getMassExcess(1000*Zl+Al);
    double Mh = Ah*amuMeV + getMassExcess(1000*Zh+Ah);
    lightFragment->mass = Ml;
    heavyFragment->mass = Mh;

    // Fission fragment kinetic energy is TKE*m_i/(m_1+m_2)
    KEl = TKE*Mh/(Ml + Mh);
    KEh = TKE*Ml/(Ml + Mh);

    // Store the sampled event
    lightFragment->A = Al;
    lightFragment->Z = Zl;
    lightFragment->N = Al-Zl;
    lightFragment->U = Ul;
    lightFragment->spin = spinl;
    lightFragment->parity = parl;
    lightFragment->KE = KEl;
    lightFragment->TKE = TKE;
    lightFragment->energyRelease = energyRelease;

    heavyFragment->A = Ah;
    heavyFragment->Z = Zh;
    heavyFragment->N = Ah-Zh;
    heavyFragment->U = Uh;
    heavyFragment->spin = spinh;
    heavyFragment->parity = parh;
    heavyFragment->KE = KEh;
    heavyFragment->TKE = TKE;
    heavyFragment->energyRelease = energyRelease;

    delete [] momentum;

    return;
}

/*******************************************************************************
 * sampleFissionFragments
 *------------------------------------------------------------------------------
 * Performs the Monte Carlo sampling of the fission fragment yields ffy to
 * select a fragment pair (Al,Zl,KEl) - (Ah,Zh,KEh) and then shares the
 * excitation energy and select the J^Pi for each fragment via the function
 * "computeFragmentInitialConditions".
 ******************************************************************************/
void FissionFragments::sampleFissionFragments(fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, Yields *ffy) {

    // Fragment properties
    double TKE;
    double maxYz;
    double KEl, KEh;
    int iA, iZ; iA = iZ = 0; // Index versions of A Z
    int j;
    int Al, Zl;
    int Ah, Zh;
    int iZmin, iZmax; iZmin = -dZ; iZmax = dZ;

    double Ul, Uh;
    double spinl, spinh;
    int parl, parh;

    double y, yz; // Yields

    double energyRelease; // Q-value for a given fragment pair

    double * momentum;
    momentum = new double [3]; // (x,y,z) components of fragment momentum

    int i = 0; // Counter for the number of (A,TKE) sampling attempts

    // Monte Carlo sampling
    do {
	i++;

	// First sample an A
	iA = int(floor(rng_cgm()*(Amax - Amin + 1) + Amin));
	y = rng_cgm()*maxYieldA; // Sample a yield

	// Check if the sampled yield is lower than the Y(A) for the sampled A
	if (y <= YA[iA]) {

	    maxYz = maxYieldZ[iA];
	    j = 0; // Counter for the number of Z sampling attempts
	    do {
		j++;
		// Sample a Z value within dZ

		iZ = int(floor(rng_cgm()*(iZmax - iZmin + 1) + iZmin));
		yz = rng_cgm()*maxYz; // Sample a yield

		// Don't sample over 10000 times for a single Z
		if (j>=10000) { cerr << "Charge sampling exceeded!" << endl; }
	    } while (yz > YZA[iZ + dZ][iA]);
	}

	// Don't sample over 10000 times for a single event
	if (i>=10000) { cerr << "[sampleFissionFragments] Sampling did not converge!" << endl; }
    } while (y > YA[iA]);

    // TODO: CAREFUL WE CAN SAMPLE AN A<ASYM AND A Z>ZSYM, SO WHICH FRAGMENT IS THIS? (P.J.)
    // -----------------------------------------------------------------------------------------------------------------------------
    // Check if we've sampled a light or heavy fragment
    if (iA<=Ac/2) {
	Al = iA;
	Zl = iZ + int(floor(Z0[iA] + 0.5));
    } else {
	Al = Ac - iA;
	Zl = Zc - iZ - int(floor(Z0[iA] + 0.5));
    }
    Zh = Zc - Zl;
    Ah = Ac - Al;
    // -----------------------------------------------------------------------------------------------------------------------------

    Ul = 0.;
    Uh = 0.;

    // Sample the TKE value
    do{
      TKE = ffy->sampleTKE(Ah); // Sample the TKE value from a Gaussian with mean <TKE>(Ah) and width s_TKE(Ah)
    	// Compute the excitation energy sharing and the spin/parity assignment
      computeFragmentInitialConditions(TKE,Zl,Al,Zh,Ah,&Ul,&Uh,&spinl,&spinh,&parl,&parh,&energyRelease,lightFragment->preMomentum,heavyFragment->preMomentum);
    } while (Ul==0. && Uh==0.); // Keep sampling until neither excitation energy is 0

    for (unsigned int i=0;i<3;i++) {
	lightFragment->postMomentum[i] = lightFragment->preMomentum[i];
	heavyFragment->postMomentum[i] = heavyFragment->preMomentum[i];
    }

    double Ml = Al*amuMeV + mass_excess(Zl,Al);
    double Mh = Ah*amuMeV + mass_excess(Zh,Ah);
    lightFragment->mass = Ml;
    heavyFragment->mass = Mh;

    // Fission fragment kinetic energy is TKE*m_i/(m_1+m_2)
    KEl = TKE*Mh/(Ml+Mh);
    KEh = TKE*Ml/(Ml+Mh);

    // Store the sampled event
    lightFragment->A = Al;
    lightFragment->Z = Zl;
    lightFragment->N = Al-Zl;
    lightFragment->U = Ul;
    lightFragment->spin = spinl;
    lightFragment->parity = parl;
    lightFragment->KE = KEl;
    lightFragment->TKE = TKE;
    lightFragment->energyRelease = energyRelease;

    heavyFragment->A = Ah;
    heavyFragment->Z = Zh;
    heavyFragment->N = Ah-Zh;
    heavyFragment->U = Uh;
    heavyFragment->spin = spinh;
    heavyFragment->parity = parh;
    heavyFragment->KE = KEh;
    heavyFragment->TKE = TKE;
    heavyFragment->energyRelease = energyRelease;

    delete [] momentum;

    return;   
}

/*******************************************************************************
 * computeFragmentInitialConditions
 *------------------------------------------------------------------------------
 * Given a total kinetic energy TKE and fragment pair (Zl,Al) - (Zh,Ah), find
 * the energy sharing and spin/parity for both fission fragments. Then solves
 * for the fragment and pre-fission neutron emission in CM and lab frames.
 ******************************************************************************/
void FissionFragments::computeFragmentInitialConditions(double TKE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh, double *spinl, double *spinh, int *parl, int *parh, double *energyRelease, double* lightFragmentMomentum, double* heavyFragmentMomentum) {

    double addToEnergyRelease, totalExcitationEnergy;

    // Compute the total excitation energy (i.e. Q = TKE + TXE so TXE = Q - TKE)
    addToEnergyRelease = incidentEnergy + SnCompound;
    *energyRelease = massExcessCompound	- getMassExcess(1000*Zl+Al) - getMassExcess(1000*Zh+Ah); // Q-value

    totalExcitationEnergy = *energyRelease + addToEnergyRelease - TKE; // TXE

    // Check for negative TXE
    if (totalExcitationEnergy<0) {

	*Ul=0.0; // Set to 0 excitation energy
	*Uh=0.0;
    } else {
	// Compute the excitation energy sharing
	computeFragmentExcitationEnergy(totalExcitationEnergy,Zl,Al,Zh,Ah,Ul,Uh,energySortingOption);

	// Compute the spins Jl and Jh
	getInitialSpin(Zl,Al,Ul,spinl);
	getInitialSpin(Zh,Ah,Uh,spinh);

	// Sample the parities
	*parl = 1; if (rng_cgm()<0.5) {*parl = -1;}
	*parh = 1; if (rng_cgm()<0.5) {*parh = -1;}
    }

    // IS: new routine that computes the fragment momentum, as well as the pre-fission 
    // neutron momenta in the CM and in the lab
    
    boost_pfn (Zc0,Ac0,incidentEnergy0,TKE,Al,Zl,Ah,Zh,Ac0-Ac,anisotropyCoefficient,preeq_flag,cmEnergy_pfn,labEnergy_pfn,neutronType_pre,cmVelocities_pfn,labVelocities_pfn,Pcostheta,Ppreeq,lightFragmentMomentum,heavyFragmentMomentum);
    
  /*  IS:  this is replaced by the above routine
    //-- fragment masses
    double Ml = Al*amuMeV+getMassExcess(1000*Zl+Al);
    double Mh = Ah*amuMeV+getMassExcess(1000*Zh+Ah);

    // compute momentum vector (px,py,pz) ----------------------------------------
    double Pf0 = sqrt( (2.0*TKE)/(1.0/Ml+1.0/Mh) );

    // choose a (random) direction for the light fragment first (isotropic emission)
    double phi = rng_cgm()*twopi;
    double costheta = 2.0*rng_cgm()-1.0;
    double sintheta = sqrt(1.0-costheta*costheta);

    pf[0] = Pf0*sintheta*cos(phi);
    pf[1] = Pf0*sintheta*sin(phi);
    pf[2] = Pf0*costheta;
*/		
    return;
}

// TODO: I DON'T THINK THIS IS USED ANYWHERE... DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************
 * sortExcitationEnergy
 *------------------------------------------------------------------------------
 * Excitation energy sorting mechanism through a ratio of temperatures RT,
 * which is provided. Takes in the total excitation energy Utot to be sorted
 * and the fragment pair Al-Ah. Assumes only the most probable Z for a given
 * fragment mass A.
 ******************************************************************************/
void FissionFragments::sortExcitationEnergy (double Utot, int Al, int Ah, double *Ul, double *Uh, double RT) {

    double alf, ahf; // Level densitiy parameters for the light and heavy fragment
    double epsilon; // Difference in alf

    double El = Utot/2.0; // First assume a simple 50-50 split
    double Eh = Utot/2.0;

    int il = 0; do {il++;} while (ldEnergyGrid[il]<El); // Find the energy index corresponding to this split
    int ih = il;
   
    // Iteratively adjust El and Eh (and il and ih) until this energy ratio matches RT = [El*ahf(Eh)]/[Eh*alf(El)]
    int count=0;
    do {

	alf = levelDensityParameters [Al][dZ][il]; // Pull the level density parameter for (Zl,Al) at energy il -- assumes Zl = Zp(Al)
	ahf = levelDensityParameters [Ah][dZ][ih]; // Same for heavy fragment

	Eh = ahf/(alf*RT*RT + ahf)*Utot; // Calculate the Eh and El corresponding to the RT value
	El = Utot - Eh;

	il = 0; do {il++;} while (ldEnergyGrid[il]<El); // Find the new light and heavy energy indices (il and ih)
	ih = 0; do {ih++;} while (ldEnergyGrid[ih]<Eh);

	epsilon = abs(alf - levelDensityParameters[Al][dZ][il]); // Check if the sample El converged towards the right alf -- works because a is monotonically increasing

	// Don't sample more than 10 times
	if (++count>10) {cerr << "[sortExcitationEnergy] did not converge!\n"; exit(-1);}

    } while (epsilon>0.01); // Cutoff in alf
   
    *Uh = Eh;
    *Ul = El;
    return;
}
// -----------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************
 * computeFragmentExcitationEnergy
 *------------------------------------------------------------------------------
 * Given a total excitation energy TXE and fragment pair (Zl,Al) - (Zh,Ah), we
 * share the excitation energy between the light and heavy fragment via the 
 * energy sorting mechanism given by "sortingOption". << more info before release >>
 ******************************************************************************/
void FissionFragments::computeFragmentExcitationEnergy(double TXE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh, string sortingOption) {

    // Intrinsic and collective excitation energy (used for sortingOption = "maxEntropy"
    double Uint;  // Total intrinsic excitation energy is shared between the fragments
    double Ucoll; // Total collective energy is the energy stored in the collective degrees of freedom at scission (not shared)

    // TODO: THIS NEEDS TO BE ADJUSTED (BACK TO Z0INIT) IF ITS EVER GOING TO BE USED (P.J.)
    int iZl = Zl;// - Z0init[Al] + dZ1; // Determine the index corresponding to Zl
    int iZh = Zh;// - Z0init[Ah] + dZ1; // and Zh

    if (sortingOption == "maxEntropy") {

	// TODO: IMPROVE THE MAX ENTROPY METHOD WITH THE DEFORMATION ENERGIES (WORK WITH P. JAFFKE, P. TALOU, J. RANDRUP, AND R. VOGT)
	// -----------------------------------------------------------------------------------------------------------------------------

	// PROCESS:
	// 1) Somewhere before, read the deformation energies V_def(A,Z,U,epsilon) -- maybe don't include U-dependence yet?
	// 2) Calculate the total potential energy at scission: V_tot = V_def(LF) + V_def(HF) + V_int(LF,HF) + V_prox(LF,HF) for a variety of epsilon_LF and epsilon_HF values
	// 3) Determine the epsilon_LF and epsilon_HF that correspond to the minimum of V_tot
	// 4) Calculate V_int = TXE - V_def(LF,epsilon_LF) - V_def(HF,epsilon_HF) and share V_int via max. entropy -- U_int(LF) + U_int(HF) = V_int
	// 5) Compute the total fragment excitation energies by U_LF = U_int(LF) + V_def(LF,epsilon_LF) and U_HF = U_int(HF) + V_def(HF,epsilon_HF)
	// X? What to do about the fact that scission deformation relaxes into g.s. deformation?

	double f;
	double Utot;

	f = 0.5; // Preliminary 50% split between collective and intrinsic energy
	Utot  = TXE - incidentEnergy - SnCompound;
	Ucoll = f*Utot;
	Uint  = (1 - f)*Utot + incidentEnergy + SnCompound;
      
	// Cf252(sf)
	if (ZAIDc == 98252) {
	    if (abs(Ah-Asym)>=25) { // For more asymmetric fragments, assume less collective energy
		Ucoll = 0.3*TXE;
	    } else {
		Ucoll = 0.6*TXE;
	    }
	}

	// U235(n,f)
	if  (ZAIDc==92236) {
	    if (abs(Ah-Asym)>=25) { // Not sure where these numbers come from?
		Ucoll = 0.2*TXE;
	    } else {
		Ucoll = 0.75*TXE;
	    }
	}

	// Pu239(n,f)
	if  (ZAIDc==94240) {
	    if (abs(Ah-Asym)>=20) {
		Ucoll = 0.7*TXE;
	    } else {
		Ucoll = 0.7*TXE;
	    }
	}
	// -----------------------------------------------------------------------------------------------------------------------------

	Uint = TXE - Ucoll;

	// Max. entropy principle states that the excitation energy is shared in the configuration that contains the maximum entropy: S = log(rho)
	// We can solve for this configuration via <U_LF> = Int_0^Uint [e * rho(Zl,Al,e) * rho(Zh,Ah,Uint-e)] de / Int_0^Uint [rho(Zl,Al,e) * rho(Zh,Ah,Uint-e)] de

	// Compute from max. entropy principle
	double numerator = 0.0; // Integral in the numerator
	double denominator = 0.0; // Integral in the denominator
	double El, Eh;
	double x;

	for (int k=1; k<=NUME; k++){

	    if (ldEnergyGrid[k]>Uint) {break;} // Integrate only up to Uint
         
	    El = ldEnergyGrid[k]; // Energy of the light fragment
	    Eh = Uint - El;

	    // Linear interpolation in S = log(rho) between the energy points
	    double yl, yh;
	    double y1, y2, x1, x2;

	    int il = 0;
	    do {
		il++;
	    } while (ldEnergyGrid[il]<El);

	    // TODO: CHANGE LEVEL DENSITY TABLES FOR ZAI INSTEAD OF A AND Z. SYNCHRONIZE THE ENERGY GRID WITH LDP? (P.J.)
	    // -----------------------------------------------------------------------------------------------------------------------------
	    y1 = log(levelDensities[Al][iZl][il - 1]);	// S(E-dE)
	    y2 = log(levelDensities[Al][iZl][il]);	// S(E)
	    x1 = ldEnergyGrid[il - 1];			// E-dE
	    x2 = ldEnergyGrid[il];			// E
	    yl = (y2 - y1)/(x2 - x1)*(El - x1) + y1;	// S(El) = S(E-dE) + [S(E) - S(E-dE)]*[El - E + dE] / [dE]
	    yl = exp(yl);				// Convert back to rho

	    // Do the same for heavy fragment
	    int ih = 0;
	    do {
		ih++;
	    } while (ldEnergyGrid[ih]<Eh);

	    y1 = log(levelDensities[Ah][iZh][ih-1]);
	    y2 = log(levelDensities[Ah][iZh][ih]);
	    x1 = ldEnergyGrid[ih-1];
	    x2 = ldEnergyGrid[ih];
	    yh = (y2-y1)/(x2-x1)*(Eh-x1)+y1;
	    yh = exp(yh);

	    // Determine the product of the level densities
	    x = yl*yh;
	    numerator += El*x; // Append to the numerator
	    denominator += x;  // and denominator
	}

	// Take the ratio to find <U_LF>
	if (denominator!=0) {
	    *Ul = numerator/denominator;
	    *Uh = Uint - *Ul;
	} else { // Sharing failed
	    *Ul = 0.0;
	    *Uh = 0.0;
	}
	return;

    } else if (sortingOption == "fixedRT" || sortingOption == "Pfactors" || sortingOption == "RTA") {
	// Everything else uses the ratio of temperatures
      
	double alf, ahf; // Level density parameters for the light and heavy fragment

	// TODO: PUT SOME EXPLANATION OF PFACTORS IN HERE - BASED ON CLOSED SHELLS BLAH BLAH BLAH
	// -----------------------------------------------------------------------------------------------------------------------------
	if (sortingOption == "Pfactors") {
	    if (Pfactors[Ah][Zh]!=0) {
		RT = 1 + 0.05*(Pfactors[Al][Zl]/Pfactors[Ah][Zh] - 1.0);
	    } else {
		RT = 1.0;
	    }
	}
	// -----------------------------------------------------------------------------------------------------------------------------

	if (sortingOption == "RTA") { // Pull the A-dependent RT value
	    RT = 1.0; // Default
	    for (unsigned int i=0;i<NUMZAID_RTA;i++) {
		if (RTAdata[i].ZAID == ZAIDt) {
		    RT = RTAdata[i].ratios[Ah];
		    break;
		}
	    }
	}

	double El = TXE/2.0; // Assume a 50-50 split
	double Eh = TXE/2.0;

	int ZAIDl = 1000*Zl + Al;
	int ZAIDh = 1000*Zh + Ah;

	alf = getLDP(ZAIDl,El);
	ahf = getLDP(ZAIDh,Eh);

	double alf0; // Last value of alf

	int count=0;
	// Iteratively adjust El and Eh (and il and ih) until this energy ratio matches RT = [El*ahf(Eh)]/[Eh*alf(El)]
	do {

	    Eh = ahf/(alf*RT*RT + ahf)*TXE; // Calculate the Eh and El corresponding to the RT value
	    El = TXE - Eh;
	    alf0 = alf;

	    alf = getLDP(ZAIDl,El);
	    ahf = getLDP(ZAIDh,Eh);

	    if (++count>10) { // Cut off the convergence at 10 (rarely fails)
		cerr << "[computeFragmentExcitationEnergy] did not converge for ZAID: " << ZAIDl << ", " << ZAIDh << endl;
//		exit(-1);
		// Defaults to 50-50 split
		*Uh = TXE/2.0;
		*Ul = TXE/2.0;
	    }

	} while (abs(alf-alf0)>0.01); // Cutoff in alf

	*Uh = Eh;
	*Ul = El;
      
    } else {
	cerr << "sorting option not implemented yet!" << endl; exit(0);
    }
    return;
}

/*******************************************************************************
 * readRTAh
 *------------------------------------------------------------------------------
 * Set tables R_T(Ah) for a given reaction. Numbers obtained from fitting CGMF
 * calculations to nubar(A) data. The data source is listed for each reaction.
 ******************************************************************************/
void FissionFragments::readRTAh(int Zt, int At) {

    int j;
    int zaidt = 1000*Zt + At;

    std::fill_n(RTAh,NUMA,1.0);

    switch (zaidt) {

	/* URANIUMS */
	case 92233: // U233
	    // from Ah=117 (sym) to Ah=159
	    static double U233_RTAh [43] = { // Adjusted to fit Nishio, 1998
		1.200, 1.200, 1.200, 1.200, 1.200, 1.200, 1.150, 1.150, 1.120, 1.120,
		1.200, 1.180, 1.180, 1.180, 1.180, 1.170, 1.110, 1.050, 1.020, 1.100,
		1.220, 1.230, 1.240, 1.240, 1.240, 1.240, 1.260, 1.280, 1.280, 1.240,
		1.240, 1.220, 1.170, 1.160, 1.160, 1.160, 1.160, 1.160, 1.160, 1.160,
		1.160, 1.160, 1.160};

	    j = -1;
	    for (int i=117; i<=159; i++) {
		RTAh[i] = U233_RTAh[++j];
	    }
	    break;

	case 92235: // U235
	    // from Ah=118 (sym) to Ah=166
	    static double U235_RTAh [49] = { // Adjusted to fit Vorobyev, 2010
		1.000, 1.096, 1.191, 1.287, 1.382, 1.478, 1.574, 1.590, 1.588, 1.388,
		1.307, 1.246, 1.195, 1.084, 1.052, 1.040, 1.047, 1.135, 1.183, 1.301,
		1.278, 1.255, 1.281, 1.228, 1.225, 1.221, 1.218, 1.214, 1.321, 1.228,
		1.204, 1.201, 1.197, 1.194, 1.191, 1.187, 1.184, 1.180, 1.177, 1.174,
		1.170, 1.167, 1.163, 1.160, 1.157, 1.153, 1.150, 1.146, 1.143 };

	    j = -1;
	    for (int i=118; i<=166; i++) {
		RTAh[i] = U235_RTAh[++j];
	    }
	    break;

	case 92238: // U238
	    // from Ah=119 (sym) to Ah=165
	    static double U238_RTAh [49] = { // Old U235(nth,f) values copied for U238(n,f)
		1.000, 1.096, 1.191, 1.287, 1.382, 1.478, 1.574, 1.669, 1.588, 1.508,
		1.427, 1.346, 1.265, 1.184, 1.192, 1.200, 1.207, 1.215, 1.223, 1.231,
		1.238, 1.235, 1.231, 1.228, 1.225, 1.221, 1.218, 1.214, 1.211, 1.208,
		1.204, 1.201, 1.197, 1.194, 1.191, 1.187, 1.184, 1.180, 1.177, 1.174,
		1.170, 1.167, 1.163, 1.160, 1.157, 1.153, 1.150, 1.146, 1.143 };
         
	    j = -1;
	    for (int i=119; i<166; i++) {
		RTAh[i] = U238_RTAh[++j];
	    }
	    break;

	/* NEPTUNIUMS */
	case 93237: // Np237
	    // from Ah=119 (sym) to Ah=156
	    static double Np237_RTAh [38] = { // Adjusted to fit Muller, 1981 (En = 0.8 MeV)
		1.200, 1.190, 1.200, 1.210, 1.200, 1.200, 1.200, 1.210, 1.190, 1.200,
		1.200, 1.200, 1.190, 1.180, 1.170, 1.170, 1.170, 1.180, 1.200, 1.200,
		1.200, 1.200, 1.210, 1.200, 1.190, 1.200, 1.200, 1.210, 1.210, 1.200,
		1.210, 1.220, 1.230, 1.230, 1.250, 1.260, 1.260, 1.260 };

	    j = -1;
	    for (int i=119; i<=156; i++) {
		RTAh[i] = Np237_RTAh[++j];
	    }
	    break;

	/* PLUTONIUMS */
	case 94239: // Pu239
	    // from Ah=120 (sym) to Ah=170
	    static double Pu239_RTAh [51] = { // Adjusted to fit Apalin, 1965
		1.000, 1.130, 1.120, 1.160, 1.140, 1.190, 1.225, 1.240, 1.330, 1.305,
		1.260, 1.110, 1.105, 1.075, 1.020, 1.010, 1.000, 1.080, 1.100, 1.075,
		1.065, 1.045, 1.030, 1.015, 1.005, 0.985, 0.955, 0.935, 0.910, 0.895,
		0.875, 0.850, 0.825, 0.850, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830,
		0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830, 0.830,
		0.830 };

	    j = -1;
	    for (int i=120; i<=170; i++) {
		RTAh[i] = Pu239_RTAh[++j];
	    }
	    break;

	case 94240: // Pu240
	    // from Ah=120(sym) to Ah=170
	    static double Pu240_RTAh[51] = { // Adjusted by P. Jaffke to match Wahl nubar(A)
		1.300, 1.250, 1.150, 1.120, 1.150, 1.200, 1.230, 1.270, 1.300, 1.310,
		1.300, 1.270, 1.200, 1.100, 1.110, 1.270, 1.290, 1.300, 1.290, 1.280,
		1.270, 1.240, 1.210, 1.210, 1.250, 1.280, 1.300, 1.320, 1.340, 1.350,
		1.380, 1.400, 1.430, 1.450, 1.480, 1.480, 1.440, 1.400, 1.380, 1.350,
		1.330, 1.300, 1.250, 1.100, 1.000, 0.980, 0.950, 0.900, 0.870, 0.850,
		0.850};

	     j = -1;
	     for (int i=120;i<=170;i++) {
		RTAh[i] = Pu240_RTAh[++j];
	     }
	     break;

	case 94241: // Pu241
	    // from Ah=121(sym) to Ah=171
	    static double Pu241_RTAh [51] = { // Adjusted by P. Jaffke to match Wahl nubar(A)
		1.05 , 1.111, 1.164, 1.212, 1.252, 1.287, 1.316, 1.339, 1.357, 1.370,
		1.378, 1.382, 1.381, 1.377, 1.369, 1.358, 1.344, 1.327, 1.307, 1.285,
		1.262, 1.237, 1.21 , 1.182, 1.154, 1.125, 1.095, 1.066, 1.037, 1.009,
		0.981, 0.955, 0.93 , 0.907, 0.885, 0.866, 0.85 , 0.836, 0.825, 0.818,
		0.815, 0.815, 0.82 , 0.829, 0.843, 0.861, 0.886, 0.916, 0.951, 0.993,
		1.041};

	    j = -1;
	    for (int i=121; i<=171; i++) {
		RTAh[i] = Pu241_RTAh[++j];
	    }
	    break;

	case 94242: // Pu242
	    // from Ah=121(sym) to Ah=171
	    static double Pu242_RTAh [51] = { // Adjusted by P. Jaffke to match Wahl nubar(A)
		1.128, 1.175, 1.217, 1.254, 1.286, 1.313, 1.335, 1.353, 1.367, 1.377,
		1.384, 1.387, 1.386, 1.383, 1.377, 1.368, 1.357, 1.343, 1.328, 1.311,
		1.292, 1.272, 1.251, 1.229, 1.206, 1.183, 1.16 , 1.136, 1.113, 1.09 ,
		1.068, 1.046, 1.026, 1.007, 0.989, 0.973, 0.959, 0.947, 0.937, 0.93 ,
		0.925, 0.924, 0.926, 0.931, 0.94 , 0.952, 0.969, 0.99 , 1.015, 1.045,
		1.08 };

	    j = -1;
	    for (int i=121; i<=171; i++) {
		RTAh[i] = Pu242_RTAh[++j];
	    }
	    break;

	/* CALIFORNIUMS */
	case 98252: // Cf252
	    // from Ah=126 (sym) to Ah=169
	    static double Cf252_RTAh [44] = { // Adjusted to fit Gook, 2014
		1.000, 1.173, 1.346, 1.519, 1.692, 1.613, 1.534, 1.455, 1.376, 1.297,
		1.307, 1.317, 1.328, 1.338, 1.348, 1.359, 1.328, 1.298, 1.267, 1.237,
		1.196, 1.166, 1.146, 1.115, 1.085, 1.054, 1.051, 0.985, 0.980, 0.953,
		0.977, 0.991, 0.996, 1.010, 1.024, 1.038, 1.029, 1.020, 1.011, 1.002,
		0.993, 0.984, 0.975, 0.496 };

	    j = -1;
	    for (int i=126; i<=169; i++) {
		RTAh[i] = Cf252_RTAh[++j];
	    }
	    break;

	case 98254: // Cf254
	    // from Ah=127 (sym) to Ah=169
	    static double Cf254_RTAh [43] = { // Same as Cf252(sf) but shifted
		1.173, 1.346, 1.519, 1.692, 1.613, 1.534, 1.455, 1.376, 1.297, 1.307,
		1.317, 1.328, 1.338, 1.348, 1.359, 1.328, 1.298, 1.267, 1.237, 1.196,
		1.166, 1.146, 1.115, 1.085, 1.054, 1.051, 0.985, 0.980, 0.953, 0.977,
		0.991, 0.996, 1.010, 1.024, 1.038, 1.029, 1.020, 1.011, 1.002, 0.993,
		0.984, 0.975, 0.496 };

	    j = -1;
	    for (int i=127; i<=169; i++) {
		RTAh[i] = Cf254_RTAh[++j];
	    }
	    break;

	/* Not supported */
	default:
	    exit(-1);
    }
    return;
}

/*********************************************************************************
 * Read R_T(Ah) from file.
 ********************************************************************************/
void FissionFragments::readRTAh (int Zt, int At, string filename) {

//    int ZAIDt;
//    ZAIDt = 1000*Zt+At;

    int Ah;
    double RT;
  
    std::fill_n (RTAh, NUMA, 1.0);
  
    ifstream dataFile;
    dataFile.open(&filename[0]);
  
    while (!dataFile.eof()) {
	dataFile >> Ah >> RT;
	RTAh[Ah] = RT;
    }

    dataFile.close();
    return;
}


/*******************************************************************************
 * readRTA
 *------------------------------------------------------------------------------
 * Reads the RT(A) tables for the supported reactions.
 ******************************************************************************/
void readRTA (void) {
	
  NUMZAID_RTA = 10;
  RTAdata = new RTAtype [NUMZAID_RTA];

  ifstream fp;
  int Amin, Amax;
  int c, n;

  string str = RTAFILE;

  // Try current directory
  fp.open(&str[0]);

  // Then system data area
  if (!fp) {
    str = datadir + str;
    fp.open(&str[0]);
  }

  if (!fp) cgmTerminateCode("R_T(A) data file not found");

  c=0;
  while(getline(fp,str)){
    if(str[0] == '#') continue;
    while (str!="") {
      RTAdata[c].ZAID = abs(atoi(str.substr(0,7).c_str()));
      Amin = atoi(str.substr(7,12).c_str());
      Amax = atoi(str.substr(12,17).c_str());
      str=str.substr(18);
      for (int i=Amin; i<=Amax; i++) {
        n=str.find(" ");
        RTAdata[c].ratios[i] = atof(str.substr(0,n).c_str());
        if (n==-1 or i==Amax) {str=""; break; }
        str=str.substr(n+1);
      }
      c++;
    }
  }

  fp.close();

  return;

}


/*******************************************************************************
 * cleanupRTA
 *------------------------------------------------------------------------------
 * Deletes the RTA dynamically allocated memory.
 ******************************************************************************/
void cleanupRTA (void) {

    delete [] RTAdata;

    return;
}


/*******************************************************************************
 * getInitialSpin
 *------------------------------------------------------------------------------
 * Given a nucleus (Z,A), an excitation energy U, and spin, sample from the spin
 * distribution: P(J) ~ (2J + 1) exp [ -(J(J + 1))/(2B^2) ]
 ******************************************************************************/
void FissionFragments::getInitialSpin(int Z, int A, double *U, double *spin) {

    int const numberSpins = 50; // Max number of spins
    double spinProba[numberSpins];
    double oddSpin;
    int j;
    double xj;
    double B2 = spinParameter2(*U, A, Z); // spin parameter B^2 (based on mom. of inertia and frag. deformation)

    oddSpin = 0.0;
    if (A%2!=0) {oddSpin = 0.5;} // For odd masses, add 1/2hbar

    double sum = 0.0;
    for (int j=0; j<numberSpins; j++) {
	xj = float(j) + oddSpin;
	spinProba[j] = (2*xj + 1)*exp(-(xj*(xj + 1))/(2*B2));
	sum += spinProba[j];
    }

    for (int j=0; j<numberSpins; j++) spinProba[j] /= sum; // Renormalize P(J)

    do  { // Sample from P(J)
	j = int(rng_cgm()*numberSpins);
    } while (rng_cgm()>spinProba[j]);

    *spin = float(j) + oddSpin; // Return the sampled spin
    return;
}

/*******************************************************************************
 * writeYieldsInFile
 *------------------------------------------------------------------------------
 * Writes the computed Y(A), Y(Z), and Y(TKE) to outputFilename
 ******************************************************************************/
void FissionFragments::writeYieldsInFile(string outputFilename) {

    ofstream outputFile;

    outputFile.open(&outputFilename[0]);

    outputFile << "# Ac = " << Ac << " ; Zc = " <<Zc << endl;

    outputFile << endl << "# Mass Distribution " << endl << endl;
    for (int i=Amin; i<=Amax; i++) {
	outputFile << i << " " << YA[i] << endl;
    }
    outputFile << endl << "# Charge Distribution " << endl << endl;
    for (int i=Zmin; i<=Zmax; i++) {
	outputFile << i << " " << YZ[i] << endl;
    }
    outputFile << endl << "# TKE Distribution " << endl << endl;
    for (int i=TKEmin; i<=TKEmax; i++) {
	outputFile << i << " " << YTKE[i] << endl;
    }

    outputFile.close();
    return;   
}

double FissionFragments::preeq_fraction(void) {

  double a0,e0,s,f0;

  int m=-1;

  for(int i=0; i < n_reactions; i++){
    if(ZAIDt==preeq_id[i]){
      m=i;   // location of the target in the preequlibrium parameters file
      break;
    }
  }

  if(m<0)
    return 0;
  else{
    a0=preeq_params[m][0];
    e0=preeq_params[m][1];
    s =preeq_params[m][2];
    f0=preeq_params[m][3];
  }

  /*
  switch(ZAIDt){

  case 93237:
    a0=11.913;
    e0=12.948;
    s=-0.011;
    f0=-0.293;
    break;
  case 92233:
     a0=11.967;
     e0=11.721;
     s=-0.012;
     f0=-0.272;

     break;

  case 92235:
    a0=11.391;
    e0=12.549;
    s=-0.010;
    f0=-0.297;
    break;

  case 92238:
    a0=8.634;
    e0=15.022;
    s=-0.007;
    f0=-0.374;
    break;

  case 94239:
    a0=11.299;
    e0=12.602;
    s=-0.010;
    f0=-0.299;
    break;

  case 94241:
    a0=10.884;
    e0=13.102;
    s=-0.009;
    f0=-0.313;
    break;
  default:
    cout << "cannot handle this reaction" << endl;
    
  }
  */

  return 1./(1+exp((a0-incidentEnergy)/e0))+s*incidentEnergy+f0;
}

/*******************************************************************************
 * buildYields
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(A,TKE) according to 'buildOption'.
 ******************************************************************************/
void FissionFragments::buildYields(void) {

    int i,j;

    SnCompound = -mass_excess(Zc,Ac) + mass_excess(Zc,Ac-1) + mass_excess(0,1);
    if (incidentEnergy == 0.0) SnCompound = 0.0; // SPONTANEOUS FISSION
    SnCompound0 = SnCompound;

    if (yieldsOption == "YATKE") {

	// Read the Y(A,TKE) file and build the Y(A) from this
	readYieldsATKE(yieldsFilename);
	buildYA();
	// Then build Y(Z|A)
    // Initialize and compute the Wahl parameters for a Z-distribution
    ZDistribution WahlParams = ZDistribution(Z0,sZ0,FZZ0,FNZ0);
    if (WahlParams.ReturnMethod() == "LegacyCGMF") {
    	WahlParams.ComputeZDistribution(Zc,Ac,incidentEnergy+SnCompound,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    } else {
    	WahlParams.ComputeZDistribution(Zc0,Ac0,incidentEnergy0+SnCompound0,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    }
	// computeWahlParametersForChargeDistribution(); // replaces buildZp(string)

	// Build the entire Y(A,Z) and corresponding Y(Z) and Y(TKE)
	buildYZA();
	buildYZ();
	buildYTKE();

    } else if (yieldsOption == "Systematics") {

	// The min and maximum A values for the fragments
	Amin = 50;
	Amax = 200;
    // Initialize and compute the Wahl parameters for a Z-distribution
    ZDistribution WahlParams = ZDistribution(Z0,sZ0,FZZ0,FNZ0);
    if (WahlParams.ReturnMethod() == "LegacyCGMF") {
    	WahlParams.ComputeZDistribution(Zc,Ac,incidentEnergy+SnCompound,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    } else {
    	WahlParams.ComputeZDistribution(Zc0,Ac0,incidentEnergy0+SnCompound0,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    }

	// computeWahlParametersForChargeDistribution(); // Calculates the midpoint Z0, width sZ0, and even-odd shell effects FNZ0, FZZ0
	// Set up discrete levels, energy bins and level densities
	ncl[0].za.setZA(Zc,Ac);
	ncl[0].max_energy = incidentEnergy + SnCompound;
	nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za); // The maximum number of neutrons that can be emitted from (Zc,Ac) with energy Q = En + Sn

	// Read the pre-fission neutron emission probability for neutron-induced fission
	if (sf_flag) {
            emissprob_max = 0.0;
	    emissprob = new double[1]; // 100% chance of 1st chance fission
	    emissprob[0] = 1.;
	} else {
	    readMultichanceFissionData();
	}

	/* FIRST, WE NEED TO KNOW IF WE EMIT A PRE-FISSION NEUTRON */

	// Initialize transmission coefficients
	try {
	    tc = new Transmission[MAX_ENERGY_BIN];
	    energy_grid1 = new double[MAX_ENERGY_BIN];
	    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
	}
	catch (bad_alloc) {
	     cgmTerminateCode("memory allocation error in excinterface");
	}
      
	// Initialize system parameters
	statSetupInitSystem(nemit,pdt);
	getRiplDiscreteLevels(ncl,nemit);

	// Read the discrete level info from RIPL for the possible daughters
	for (int i=0 ; i<=nemit ; i++) {
	    statFixDiscreteLevels(&ncl[i]);
	    statSetupEnergyBin(&ncl[i]);
	    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp); // Get the level density parameters for this nucleus
	    statClearPopulation(&ncl[i]);
	}

	// first determine whether we have sampled a pre-equilibrium neutron
	// then calculate spectra
    	do {
	   n_pfn = getNumberPrefissionNeutrons(); // sample the number of pre-fission neutrons
    	} while (n_pfn<0);
    	if (n_pfn==0) preeq_flag=false; // catch for not sampling any pre-fission neutrons

	if (n_pfn > 0){

	// Create the pre-fission neutron spectrum sp_pfn
	sp1    = new double[MAX_ENERGY_BIN];
	sp     = new double[MAX_ENERGY_BIN];

	// Fraction of pre-equilibrium spectrum

	double f_pe=preeq_fraction();
	// initial parameterization
	//double f_pe = 1.0/(1.0+ exp((12.4924628980014-incidentEnergy)/10.2096843664751)) -0.00422967345976505*incidentEnergy-0.249024419940462;
	if (f_pe<0.) f_pe = 0.;

    // determine whether we've sampled a pre-equilibrium neutron or not 
	double rand = rng_cgm();

	//preeq_flag = false;
	if (rand<=f_pe){
	   // set flag for pre_eq neutron
	   preeq_flag = true;
	} else{
	   preeq_flag = false;
	}

	if (preeq_flag){
	   // construct the evaporation spectrum - because we need tc
	   evapInterface(pdt+neutron,tc,0,0,sp1); // Create the evaporation spectrum

	   // construct the pre-equilibrium energy spectrum
	   double *spc[MAX_CHANNEL];
	   for (int j=0; j<MAX_CHANNEL; j++) {
	      spc[j] = new double[MAX_ENERGY_BIN];
	   }
	   excitonInterface(Zc, Ac-1, incidentEnergy,spc); // Calculate the pre-equilibrium spectrum
	   for (int i=0; i<MAX_ENERGY_BIN; i++){
	      sp1[i] = spc[neutron][i];  // now we can just sample from sp1 for pre-eq or compound
     	   }
	    for (int j=0; j<MAX_CHANNEL; j++) {
		delete [] spc[j];  // this is no longer needed
	    }
	} else {
	   // construct the evaporation spectrum
	   evapInterface(pdt+neutron,tc,0,0,sp1); // Create the evaporation spectrum
	}

	// Create the energy grid
	for (int i=0;i<MAX_ENERGY_BIN;i++) {
	    energy_grid1[i] = tc[i].ecms;
	}

	}

	// TODO: IS THIS A REPEAT OF ABOVE? (P.J.)
	// -----------------------------------------------------------------------------------------------------------------------------
	ncl[0].za.setZA(Zc,Ac);
	ncl[0].max_energy = incidentEnergy + SnCompound;

	// Initialize system parameters */
	statSetupInitSystem(nemit,pdt);
	getRiplDiscreteLevels(ncl,nemit);

	for (int i=0 ; i<=nemit ; i++) {
	    statFixDiscreteLevels(&ncl[i]);
	    statSetupEnergyBin(&ncl[i]);
	    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
	    statClearPopulation(&ncl[i]);
	}
	// -----------------------------------------------------------------------------------------------------------------------------

    } else if (yieldsOption == "YAZTKEfile") {
	readYieldsAZTKE(yieldsFilename);
	buildYA();
	buildYZ();
	buildYTKE();

    } else { // Right now we only have Systematics and YATKE

	cerr << "Fission fragment build option not recognized!" << endl;
	exit(-1);
    }

    // If we built the yields, they are already symmetric
    if (yieldsOption == "Systematics") return;

    // Else we need to ensure they are symmetric
    for (i=Amin; i<Asym; i++) {
	YA[Ac-i] = YA[i];
	for (j=Zmin; j<=Zmax; j++) {
	    YZA2[Zc-j][Ac-i] = YZA2[j][i];
	}
	for (j=0; j<=2*dZ; j++) {
	    YZA[2*dZ-j][Ac-i] = YZA[j][i];
	}
    }

    // Find maximum Y(Z|A) for each A -- for sampling
    for (i=0; i<NUMA; i++) {
	maxYieldZ[i] = 0.0;
	for (j=0; j<NUMZ; j++) {
	    if (YZA2[j][i]>maxYieldZ[i]) {maxYieldZ[i] = YZA2[j][i];}
	}
    }

    // Find maximum Y(A,TKE) -- for sampling
    maxYieldATKE=0.0;
    for (int i=Amin; i<=Amax; i++) {
	for (int j=TKEmin; j<=TKEmax; j++) {
	    if (YATKE[i][j]>maxYieldATKE) {maxYieldATKE = YATKE[i][j];}
	}
    }
    return;   
}

/*******************************************************************************
 * buildYA
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(A) from Y(A,TKE)
 ******************************************************************************/
void FissionFragments::buildYA(void) {

    float sum;
    sum = 0.0;
    for (int i=0; i<NUMA; i++) {
	YA[i] = 0.0;
	for (int j=0; j<NUMTKE; j++) {
	    YA[i] = YA[i] + YATKE[i][j]; // Sum across TKE
	}
	sum = sum + YA[i];
    }
    if (sum!=0) {
	for (int i=0; i<NUMA; i++) {
	    YA[i] /= sum; // Renormalize Y(A)
	}
    }
    return;   
}

/*******************************************************************************
 * buildYA
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(A,TKE) given a yields parameterization
 * ffy
 ******************************************************************************/
void FissionFragments::buildYATKE(Yields * ffy) {

    // Initialize the min and max TKE -- for sampling
    TKEmax=-1; TKEmin=500;
    for (int a=Amin; a<Amax; ++a) {
	for (int t=0; t<NUMTKE; ++t) {
	    double y = ffy->yieldATKE(a, (double) t); // Returns the Y(A,TKE)
	    if (y>1.e-9) { // A yield cutoff
		YATKE[a][t] = y;
		if (t<TKEmin) {
		    TKEmin = t; // Reset minTKE
		}
		if (t>TKEmax) {
		    TKEmax = t; // Reset maxTKE
		}
	    } else {
		YATKE[a][t] = 0.; // Set Y(A,TKE) to 0 for Y(A,TKE) < 1e-9
	    }
	}
    }
    return;
}

/*******************************************************************************
 * buildYZA
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(Z,A) where the charge distribution uses
 * the Wahl systematics after the parameters have been set via the function
 * computeWahlParametersForChargeDistribution() -- See "Systematics of Fission
 * Product Yields" - LA-13928 by Wahl (2002)
 ******************************************************************************/
void FissionFragments::buildYZA(void) {

    double c, Zp0;
    int iZp0, iZ;
    int i, j;
    double oddeven, V, W;
   
    for (i=Amin; i<=Amax; i++) {
	c    = sZ0[i]; // Set the width c, mean Zp, and the integer corresponding to Zp
	Zp0  = Z0[i];
	iZp0 = (int) floor(Zp0+0.5);

	oddeven = 1.0; // Default odd-even effect
	for (j=0; j<=2*dZ; j++) {

	    // Loop through the Z values for a given A (iZ is the true Z value and j is the index)
	    iZ = iZp0+j-dZ;
	    YZA[j][i]   = 0.0; // Contains the Y(Z,A) for a charge index j and mass i
	    YZA2[iZ][i] = 0.0; // Contains the Y(Z,A) for a charge iZ and mass i

	    if (c>0.0) {

		// Calculate the total odd-even effect F_oe
		if (iZ%2==0) {
		    if (i%2==0) {
			oddeven = FZZ0[i]*FNZ0[i]; // F_oe = F_Z*F_N (even-even)
		    } else {
			oddeven = FZZ0[i]/FNZ0[i]; // F_oe = F_Z/F_N (even-odd)
		    }
		} else {
		    if ((i+1)%2==0) {
			oddeven = FNZ0[i]/FZZ0[i]; // F_oe = F_N/F_Z (odd-even)
		    } else {
			oddeven = 1.0/(FZZ0[i]*FNZ0[i]); // F_oe = 1/F_Z/F_N (odd-odd)
		    }
		}

		V = (iZ-Zp0+0.5)/(sqrt(2.0)*c);
		W = (iZ-Zp0-0.5)/(sqrt(2.0)*c);

		YZA[j][i]   = oddeven * 0.5 * ( erf(V) - erf(W) );
		YZA2[iZ][i] = YZA[j][i];

	    } else { // If the width is negative, just do the mean Zp
		if (j==dZ) {
		    YZA[j][i]   = 1.0;
		    YZA2[iZ][i] = 1.0;
		}
	    }
	} /* END Z LOOP */

	// Now renormalize and form Y(Z,A) = Y(A) * Y(Z|A)
	double x=0;
	for (j=0; j<=2*dZ; j++) {x+=YZA[j][i];} // Find Sum of Y(Z|A) over Z -- NOT NECESSARILY 1 BECAUSE OF FINITE dZ RANGE
	if (x!=0) {
	    for (j=0; j<=2*dZ; j++) {
		YZA[j][i] = YZA[j][i]*YA[i]/x; // Calculate Y(Z,A) for a charge index j and mass i
	    }
	}

	// Do the same for the Y(Z,A) for charge j and mass i
	double y=0;
	for (j=0; j<NUMZ; j++) {y+=YZA2[j][i];}
	if (y!=0) {
	    for (j=0; j<NUMZ; j++) {
		YZA2[j][i] = YZA2[j][i]*YA[i]/y;
	    }
	}
      
    } /* END A LOOP */
    return;
}


/*******************************************************************************
 * buildYZ
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(Z) from the Y(Z,A) distribution via
 * summing over A values for a given Z. Assumes buildYZA has been called.
 ******************************************************************************/
void FissionFragments::buildYZ(void) {

    int i,j;
    double sum = 0.0;

    // Determine min and max values of Zp(A) over all A
    Zpmin = 100;
    Zpmax = 0;

    for (int i=Amin; i<=Amax; i++) {
	if (YA[i]!=0.0 && Z0[i]<Zpmin) {Zpmin = (int)floor(Z0[i] + 0.5);}
	if (YA[i]!=0.0 && Z0[i]>Zpmax) {Zpmax = (int)floor(Z0[i] + 0.5);}
    }

    // Validate Zpmin and Zpmax
    if (Zmin>Zpmin-dZ) {cout<<"ERROR [buildYZ] Zmin > Zp(Amin)-dZ"; exit(-1);}
    if (Zmax<Zpmax+dZ) {cout<<"ERROR [buildYZ] Zmax < Zp(Amax)+dZ"; exit(-1);}

    // Loop over the Zmin to Zmax (determined from computeWahlParametersForChargeDistribution)
    for (i=Zmin; i<=Zmax; i++) {
	YZ[i] = 0.0;
	sum   = 0.0;
	for (j=Amin; j<=Amax; j++) {
	    YZ[i] += YZA2[i][j]; // Sum over A values
	}
	sum += YZ[i];
    }
    if (sum!=0) { // Renormalize Y(Z) to 1
	for (i=Zmin; i<=Zmax; i++) {
	    YZ[i] /= sum;
	}
    }
    return;   
}

/*******************************************************************************
 * buildYTKE
 *------------------------------------------------------------------------------
 * Builds the fission fragment yields Y(TKE) from the Y(A,TKE) distribution via
 * summing over A values for a given TKE. Assumes the Y(A,TKE) has been built.
 ******************************************************************************/
void FissionFragments::buildYTKE(void) {

    int i,j;
    double sum = 0.0;

    for (i=TKEmin; i<=TKEmax; i++) {
	YTKE[i] = 0.0;
	for (j=Amin; j<=Amax; j++) {
	    YTKE[i] += YATKE[j][i]; // Sum over A values
	}
	sum += YTKE[i];
    }
    for (i=TKEmin; i<=TKEmax; i++) { // Renormalize Y(TKE) to 1
	YTKE[i] /= sum;
    }
    return;
}

// TODO: I DON'T THINK THIS IS USED ANYWHERE. DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * buildYATKEfromSystematics
 *------------------------------------------------------------------------------
 * Build the yields Y[A][TKE] using systematics for <TKE>(A) and <sigma_TKE>(A).
 * It follows the form: Y(A,TKE) = Y(A)xP(TKE|A) where P(TKE) is a Gaussian
 * function whose parameters are mass-dependent.
 ******************************************************************************/
void FissionFragments::buildYATKEfromSystematics (void) {
  
  std::fill_n(meanTKE,NUMA,187.0); // default is constant
  std::fill_n(sigTKE,NUMA,10.0);   // default is constant
  
  buildSystematicsTKEvsMass (); // Zeke Blaine, 4/9/2012
  
  for (int i=Amin; i<=Amax; i++) {
    if (YA[i]!=0.0) {
      for (int j=TKEmin; j<=TKEmax; j++) {
        YATKE[i][j]= YA[i]*exp(-pow(j-meanTKE[i],2)/(2*pow(sigTKE[i],2)));
      }
    }
  }
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: I DON'T THINK THIS IS USED ANYHWERE. DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * buildYATKEfromSystematics << 1.0.6 >>
 *------------------------------------------------------------------------------
 * Build the yields Y[A][TKE] using systematics for <TKE>(A) and <sigma_TKE>(A).
 * It follows the form: Y(A,TKE) = Y(A)xP(TKE|A) where P(TKE) is a Gaussian
 * function whose parameters are mass-dependent.
 * I. Stetcu: this systematics use data from thermal, scalled to the correct
 * incident energy
 ******************************************************************************/
void FissionFragments::buildYATKEfromSystematics (Yields *ffy) {
  
  TKEmin=999; TKEmax=0;
  maxYieldATKE=0.0;
  double s=0.;
  double yieldTKE[NUMTKE];
  for (int i=Amin; i<=Amax; i++) {
    if (YA[i]>0.0) {
      std::fill_n(yieldTKE,NUMTKE,0.); // default is constant
      for (int j=0; j<NUMTKE; j++) {
        yieldTKE[j]= ffy->yieldTKE(i, j );
        if( yieldTKE[j]>0.0){
          if(TKEmin>j) TKEmin=j ;
          if(TKEmax<j) TKEmax=j ;
        }
      }
      double ss=0.;
      for(int j=0; j<NUMTKE;j++)
        ss+=yieldTKE[j];
      for(int j=0; j<NUMTKE;j++){
        YATKE[i][j]=YA[i]*yieldTKE[j]/ss;
        s+=YATKE[i][j];
      }
    }
  }
  
  for(int i=Amin; i<=Amax; i++){
    for(int j=TKEmin; j<=TKEmax;j++){
      YATKE[i][j]/=s;
      if( maxYieldATKE<YATKE[i][j] ) maxYieldATKE = YATKE[i][j];
    }
  }
  
  cout << "maxYieldATKE=" << maxYieldATKE << endl;
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: I DON'T THINK THIS IS USED ANYHWERE. DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * buildWahlYA
 *------------------------------------------------------------------------------
 * Builds the mass yields Y(A) from Wahl's systematics.
 *------------------------------------------------------------------------------
 * Reference:
 * A.C.Wahl, "Systematics of Fission-Product Yields", Los Alamos Technical
 * Report, LA-13928 (2002).
 *------------------------------------------------------------------------------
 * From Zeke Blain, blaine2@rpi.edu, F95 original version (03/28/2012).
 ******************************************************************************/
void FissionFragments::buildWahlYA (void) {
  
  // Initialize Variables
  // double del7[4], Y1, Y3, Y5, D[7], del2, del4; // Ask Zeke
  double T, PE, NTtot, A4tot, Abar;
  double Y[7], sigma[7], delta[7];
  double sig1P[3], del5P[3], NTP[3];
  
  double Yield2[NUMA];
  
  int stop;
  int D1, F, G, H, J, K, M;
  double ECOR, NT1, Ip, Cp, Sp, Id, Cd, Sc, Dh, Sh, SYM, TH, Fz, Fn, nus1, Rf, Rj, sum, sum2;
  
  double Nus[NUMA], Nuh[NUMA], NuL[NUMA], R[NUMA], Nut[NUMA];
  
  // allocate (yields%Ya(Amin:Amax))
  // yields%Ya = 0.0_rk
  
  // Input parameters for calculating FN,FZ,sigz,delz
  double Sn = -getMassExcess(1000*Zc+Ac)+getMassExcess(Zc*1000+Ac-1)+getMassExcess(1);
	PE = incidentEnergy + Sn;
  
  double sig1[6] = {2.808, 8.685, -0.0454, 0.372, -0.620, -0.0122};
  double sig6[4] = {3.17, 0.0, 0.303, 0.0};
  double del5[6] = {25.34, 18.55, -0.0402, -1.220, -1.732, 0.0};
  double A4[4]   = {136.66, -0.177, 0.060, -0.038};
  double Y2[4]   = {43.00, -1.91, -3.41, 0.0};
  double NT[6]   = {1.563, 16.66, -0.00804, 0.0918, 0.0, 0.0};
  
  // Calculate Values for nu
  double En = incidentEnergy;
  Ip=20.0;
  if (Zt==98 && Ac==252) Ip=0;
  Cp=72.17;
  Id=-5.0;
  Cd=18.17;
  //     D1=78
  //     F=104
  //     G=107
  //     H=117
  //     J=127
  //     K=130
  //     M=156
  Fz=pow(-1.0,Zt);
  Fn=pow(-1.0, Ac-Zt);
  TH=11.47-0.166*Zt*Zt/Ac+0.093*(2-Fz-Fn)-Sn;
  if (Zt==98 && Ac==252) {
    NT1=2.286+0.147*(Zt-92)+0.054*(Ac-236)+0.040*(2-Fz-Fn)+(0.145-0.0043*(Ac-236));
  } else {
    NT1=2.286+0.147*(Zt-92)+0.054*(Ac-236)+0.040*(2-Fz-Fn)+(0.145-0.0043*(Ac-236))*(En-TH);
  }
  ECOR=exp(-0.05*PE);
  if ((NT1+1)>4) {
    nus1=NT1+1;
  } else {
    nus1=4.0;
  }
  Sp=(Ac-nus1)/2;
  Dh=Id/sqrt(3.14159*Cd);
  SYM=(Ac-nus1)/2.0;
  Sh=Id*(exp(-1*pow(130.0-SYM,2)/Cd))/sqrt(3.14159*Cd);
  Sc=Ac-130.0-(NT1-Dh-Sh);
  Rf=0.9-0.0075*PE;
  Rj=0.1+0.0075*PE;
  for (int i=Amin; i<=Amax; i++) {
    Nus[i] = Ip*exp(-1*pow(i-Sp,2)/Cp)/sqrt(3.14159*Cp);
    Nuh[i] = Id*exp(-1*pow(i-130.,2)/Cd)/sqrt(3.14159*Cd);
    NuL[i] = Id*exp(-1*pow(i-Sc,2)/Cd)/sqrt(3.14159*Cd);
    Nut[i] = (NT1+ECOR*(Nus[i]+Nuh[i]+NuL[i]));
  }
  D1 = (int) (Ac-156-Nut[156]);
  F = (int) (Ac-135-Nut[135]);
  G = (int) (Ac-130-Nut[130]);
  H = (int) ((128+(Ac-128-Nut[128]))/2);
  J = 130;
  K = 135;
  M = 156;
  for (int i=Amin; i<=Amax; i++) {
    if (i<D1) {
      R[i]=0.2;
    } else if (i<=F && i>=D1) {
      R[i]=0.2+(i-D1)*(Rf-0.2)/(F-D1);
    } else if (i<G && i>F) {
      R[i]=Rf;
    } else if (i<=H && i>=G) {
      R[i]=0.5+(H-i)*(Rf-0.5)/(H-G);
    } else if (i<J && i>H) {
      R[i]=0.5-(i-H)*(0.5-Rj)/(J-H);
    } else if (i<=K && i>=J) {
      R[i]=Rj;
    } else if (i<M && i>K) {
      R[i]=Rj+(i-K)*(0.8-Rj)/(M-K);
    } else {
      R[i]=0.8;
    }
    Nut[i] = R[i]*(NT1+ECOR*(Nus[i]+Nuh[i]+NuL[i]));
  }
  
  // Calculate input parameter for any given fissioning system
  for (int i=0; i<3; i++) {
    sig1P[i] = sig1[i]+sig1[i+3]*(Zt-92);
    del5P[i] = del5[i]+del5[i+3]*(Zt-92);
    NTP[i]   = NT[i]+NT[i+3]*(Zt-92);
  }
  sigma[0] = sig1P[0]+(sig1P[1]-sig1P[0])*(1-exp(sig1P[2]*PE));
  sigma[1] = 2.45;
  sigma[2] = 8.6;
  if (sigma[0]>sigma[2]) sigma[2]=sigma[0];
  sigma[3] = sigma[1];
  sigma[4] = sigma[0];
  sigma[5] = sig6[0]+sig6[2]*(Zt-92);
  sigma[6] = sigma[5];
  NTtot = NTP[0]+(NTP[1]-NTP[0])*(1-exp(NTP[2]*PE));
  A4tot = (A4[0]+A4[2]*(Zt-92)+(A4[1]+A4[3]*(Zt-92))*PE);
  A4tot = A4tot-Nut[int(A4tot)];
  delta[4] = del5P[0]+(del5P[1]-del5P[0])*(1-exp(del5P[2]*PE));
  delta[0] = -delta[4];
  delta[3] = A4tot-(Ac-NT1)/2;
  delta[1] = -delta[3];
  delta[2] = 0;
  delta[6] = 30.31;
  delta[5] = -delta[6];
  Abar = (Ac-NT1)/2;
  if (PE>8.0 && PE<=20.0) {
    Y[5] = 6.8-(6.8/12.0)*(PE-8.0);
  } else if (PE<8.0) {
    Y[5] = 6.8;
  } else if (Zt==93) {
    Y[5] = Y[5]/2.0;
  } else if (Zt<93 || PE>20.0) {
    Y[5] = 0;
  }
  Y[6] = Y[5];
  Y[1] = Y2[0]+Y2[2]*(Zt-92)+(Y2[1]+Y2[3]*(Zt-92))*PE;
  if (Y[1]<0) Y[1]=0;
  Y[3] = Y[1];
  if (PE<11.96) {
    Y[2] = 4.060*exp(0.470*(PE-11.96));
  } else {
    T = -0.030+0.0050*(Ac-236);
    Y[2] = 4.060+86.02*(1.0-exp(T*(PE-11.96)));
  }
  Y[0]=0.0;
  Y[4]=Y[0];
  stop=0;
  sum2=0;
  while (stop==0) {
    sum=0;
    for (int i=Amin; i<=Amax; i++) {
      Yield2[i]=0;
      for (int j=0; j<7; j++) {
        Yield2[i] = Yield2[i]+Y[j]*exp(-(pow(i-Nut[i]-Abar+delta[j],2))/(2*pow(sigma[j],2)))/(sigma[j]*sqrt(2*3.14159));
      }
    }
    for (int i=Amin; i<Amax; i++) {
      sum = sum+min(Yield2[i],Yield2[i+1])+abs(Yield2[i]-Yield2[i+1])/2.0;
    }
    if (sum>=199.5 && sum<=200.5) {
      stop=1;
    } else {
      Y[0]=Y[0]+0.1;
      Y[4]=Y[0];
    }
  }
  for (int i=Amin; i<=Amax; i++) {
    Yield2[i]=0;
    for (int j=0; j<7; j++) {
      Yield2[i] = Yield2[i]+Y[j]*exp(-(pow(i-Nut[i]-Abar+delta[j],2))/(2*pow(sigma[j],2)))/(sigma[j]*sqrt(2*3.14159));
    }
  }
  
  for (int i=Amin; i<=Amax; i++) {
    YA[i] = Yield2[i];
    sum2 += Yield2[i];
  }
  for (int i=Amin; i<=Amax; i++) {
    YA[i] /= sum2;
  }
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: I DON'T THINK THIS IS USED ANYHWERE. DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 *                   SUBROUTINE buildSystematicsTKEvsMass
 *------------------------------------------------------------------------------
 * based on work by Ezekiel Blaine, RPI, blaine2@rpi.edu, 4/9/2012.
 * Using Hambsch data for U-235, Pu-239 and Cf-252.
 ******************************************************************************/
void FissionFragments::buildSystematicsTKEvsMass (void) {
  
  int i;
  double A1, A2, A3, A4, A5, P1;
  
  P1=1.00+0.0112*(Zc-94)+0.0009*(Ac-240)-0.00135*pow(float(Zc-94),2);
  A3=15.20048736;
  A4=-0.0798617477;
  A5=0.0001558700726;
  A1=39639.39298*(1+0.000045*(Zc-94)+0.000022*(Ac-240)-0.0000018*pow(float(Ac-240),2));
  A2=-1272.328785*(1+0.0000013*(Zc-94)+0.0000006*(Ac-240)-0.000000055*pow(float(Ac-240),2));
  
  meanTKE[Asym]=A1+A2*Asym+A3*pow(float(Asym),2)+A4*pow(float(Asym),3)+A5*pow(float(Asym),4);
  
  for (i=1; i<=15; i++) {
    meanTKE[Asym+i] = A1+A2*(Asym+i)+A3*pow(float(Asym+i),2)+A4*pow(float(Asym+i),3)+A5*pow(float(Asym+i),4);
    meanTKE[Asym-i] = meanTKE[Asym+i];
  }
  
  for (i=21; i<=Amax-Asym; i++) {
    meanTKE[Asym+i] = 324.5*P1-1.028*(Asym+i);
    meanTKE[Asym-i] = meanTKE[Asym+i];
    if (Asym-i<0) { cerr << "[buildSystematicsTKEvsMass] error; negative array index!\n"; exit(-1); }
  }
  
  for (i=16; i<=20; i++) {
    meanTKE[Asym+i] = (meanTKE[Asym+21]-meanTKE[Asym+15])/(21-15)*(i-15)+meanTKE[Asym+15];
    meanTKE[Asym-i] = meanTKE[Asym+i];
  }
  
}
// -----------------------------------------------------------------------------------------------------------------------------

// TODO: I DON'T THINK THIS IS USED ANYHWERE. DELETE? (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * checkDistributions
 *------------------------------------------------------------------------------
 * Check the distributions of fission fragments produced by the function
 * generateInitialFissionFragmentHistories.
 ******************************************************************************/
void FissionFragments::checkDistributions(string inputFilename, string outputFilename) {

    ifstream inputFile;
    ofstream outputFile;

    // TODO: WHY IS THERE A DIFFERENT VARIABLE HERE numberEnergies? SHOULDN'T IT BE THE SAME AS WHAT'S IN CONFIG-FF.H? (P.J.)
    // -----------------------------------------------------------------------------------------------------------------------------
    //  const int NUME = 501; // number of points on energy grid
    const int numberEnergies = 401; // number of points on energy grid
    // -----------------------------------------------------------------------------------------------------------------------------

    int i;

    int iTKE, Al, Ah, Zl, Zh;
    double spinl,spinh;

    double Ul, Uh;
    double KEl, KEh, TKE;
    double Q;

    double massYields[NUMA];
    double chargeYields[NUMZ];
    double tkeYields[NUMTKE];

    double excitationEnergies[NUMA];
    double spins[NUMA];

    double sTKEvsMass[NUMA];
    double TKEvsMass[NUMA];
    double TXEvsMass[NUMA];

    double initialEnergiesLF[numberEnergies];
    double initialEnergiesHF[numberEnergies];
    double initialTXE[numberEnergies];

    // Initialize
    double **YieldsATKE;
    YieldsATKE = new double *[NUMA];
    for (int i=0; i<NUMA; i++) {
	YieldsATKE[i] = new double [NUMTKE];
	for (int j=0; j<NUMTKE; j++) {
	    YieldsATKE[i][j] = 0.0;
	}
    }

    int numberEventsPerMass[NUMA];
    int numberEventsPerCharge[NUMZ];
    int numberEventsPerTKE[NUMTKE];

    int totalNumberEvents;
    bool isSymmetric = false;

    for (int i=0; i<NUMA; i++) {
	excitationEnergies[i]  = 0.0;
	spins[i]               = 0.0;
	TKEvsMass[i]           = 0.0;
	TXEvsMass[i]           = 0.0;
	numberEventsPerMass[i] = 0;
	massYields[i]          = 0.0;
    }

    std::fill_n(chargeYields,NUMZ,0);
    std::fill_n(tkeYields,NUMTKE,0);

    std::fill_n(numberEventsPerCharge,NUMZ,0);
    std::fill_n(numberEventsPerTKE,NUMTKE,0);

    std::fill_n(initialEnergiesLF,numberEnergies,0.0);
    std::fill_n(initialEnergiesHF,numberEnergies,0.0);
    std::fill_n(initialTXE,numberEnergies,0.0);
    std::fill_n(sTKEvsMass,numberEnergies,0.0);

    double **ZAYields;
    ZAYields = new double * [NUMZ];
    for (int i=0; i<NUMZ; i++) {
	ZAYields[i] = new double [NUMA];
	std::fill_n (ZAYields[i], NUMA, 0.0);
    }

    // TODO: AGAIN, SHOULDN'T THIS BE DEFINED IN CONFIG-FF.H? (P.J.)
    // -----------------------------------------------------------------------------------------------------------------------------
    // define energy grid
    double energyGrid[numberEnergies];
    double deltaE = 0.25;
    // -----------------------------------------------------------------------------------------------------------------------------
    for (i=0; i<numberEnergies; i++) {
	energyGrid[i] = i*deltaE;
    }

    totalNumberEvents = 0;
    inputFilename = datadir + inputFilename;
    inputFile.open(&inputFilename[0]); // Load the history file
    if (!inputFile) { cerr << "[checkDistributions] ERROR: cannot find history file!" << endl; exit(-1); }
    double avTKE = 0.;

    Amin   = 1000;
    Amax   = 0;
    TKEmin = 1000;
    TKEmax = -1;
    // Loop through and read the history entries and add them to the distributions

    int Z, A, parity;
    double KE, U, spin;

    while (!inputFile.eof()) {

	inputFile >> Z >> A >> KE >> U >> spin >> parity;

	if (A<Asym) {
	    Ul = U;
	} else {
	    Uh = U;
	}

	excitationEnergies[A] += U;
	spins[A] += spin;

	// Verify that the energies aren't outside of the grid
	i = int(U/deltaE) + 1;
	if (i>numberEnergies) {
	    cout << "U = " << U << "\n";
	    cerr << "[checkDistributions] ERROR: should increase numberEnergies!" << endl;
	    exit(-1);
	}

	if (A<Asym) {
	    initialEnergiesLF[i]++;
	} else {
	    initialEnergiesHF[i]++;
	}

	if (A>=Asym) {
	    i = int((Ul+Uh)/deltaE)+1;
	    initialTXE[i]++;
	    TXEvsMass[A] += (Ul+Uh);
	}

	numberEventsPerMass[A]++;
	numberEventsPerCharge[Z]++;

	massYields[A]++;
	chargeYields[Z]++;
	totalNumberEvents++;

	ZAYields[Z][A]++;

	if (A>=Asym && !isSymmetric) {
	    KEl = KE;
	    Al = A;
	    Zl = Z;
	    if (A==Asym) isSymmetric = true;
	} else {
	    KEh = KE;
	    Ah = A;
	    Zh = Z;
	    TKE = KEl + KEh;
	    iTKE = (int)floor(TKE+0.5);
	    numberEventsPerTKE[iTKE]++;

	    tkeYields[iTKE]++;
	    TKEvsMass[Al]+= TKE;
	    TKEvsMass[Ah]+= TKE;

	    YieldsATKE[Al][iTKE]++;
	    YieldsATKE[Ah][iTKE]++;

	    // Q-values
	    Q = (mass_excess(Zc,Ac) - mass_excess(Zl,Al) - mass_excess(Zh,Ah) + incidentEnergy+SnCompound);

	    QfA[Al] += Q;
	    QfA[Ah] += Q;
	    QfTKE[iTKE] += Q;

	    if (A==Asym) {isSymmetric = false;} // Reset
	}

	if (A>=Asym) {Ul = 0.0;} else {Uh = 0.0;} // Reset
    }
    inputFile.close();
    // Two lines = one event
    totalNumberEvents /= 2;

    // Printing ------------------------------------------------------------------------------
    outputFilename = datadir + outputFilename;
    outputFile.open(&outputFilename[0]);

    int gnuplotIndex = -1;

    // Mass Yields Y(A)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] Mass Yields Y(A)" << endl << endl;
    for (int i=0; i<NUMA; i++) {
	if (numberEventsPerMass[i]!=0) {
	    massYields[i] /= totalNumberEvents;
	    outputFile << i << " " << massYields[i] << endl;
	}
    }
    outputFile << endl;

    // Charge Yields Y(Z)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] Charge Yields Y(Z)" << endl << endl;
    for (int i=0; i<NUMZ; i++) {
	if (numberEventsPerCharge[i]!=0) {
	    chargeYields[i] /= totalNumberEvents;
	    outputFile << i << " " << chargeYields[i] << endl;
	}
    }
    outputFile << endl;

    // TKE Yields Y(TKE)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] TKE   Y(TKE)  <Qf>(TKE)" << endl << endl;
    for (int i=0; i<NUMTKE; i++) {
	if (numberEventsPerTKE[i]!=0) {
	    tkeYields[i] /= totalNumberEvents;
	    YTKEf[i] /= totalNumberEvents;
	    QfTKE[i] /= numberEventsPerTKE[i];
	    outputFile << i << " " << tkeYields[i] << " " << QfTKE[i] << " " << YTKEf[i] << endl;
	}
    }
    outputFile << endl;

    // P(U)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] Initial Excitation Energy Probability P(Ui)" << endl << endl;
    double sum = 0, sumLF = 0, sumHF = 0;
    for (int i=0; i<numberEnergies; i++) {
	sum   += initialTXE[i];
	sumLF += initialEnergiesLF[i];
	sumHF += initialEnergiesHF[i];
    }
    sum   *= deltaE;
    sumLF *= deltaE;
    sumHF *= deltaE;
    for (int i=0; i<numberEnergies; i++) {
	initialTXE[i]        /= sum;
	initialEnergiesLF[i] /= sumLF;
	initialEnergiesHF[i] /= sumHF;
	outputFile << energyGrid[i] << " " << initialEnergiesLF[i] << " " <<
	initialEnergiesHF[i] << " " << initialTXE[i] << endl;
    }
    outputFile << endl;

    // <U>=f(A)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] <U>=f(A) " << endl << endl;
    for (int i=0; i<NUMA; i++) {
	if (numberEventsPerMass[i]!=0) { excitationEnergies[i] /= numberEventsPerMass[i]; }
	    outputFile << i << " " << excitationEnergies[i] << endl;
    }
    outputFile << endl;

    // <spin>=f(A)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] <spin>=f(A) " << endl << endl;
    for (int i=0; i<NUMA; i++) {
	if (numberEventsPerMass[i]!=0) {
	    spins[i] /= numberEventsPerMass[i];
	    outputFile << i << " " << spins[i] << endl;
	}
    }
    outputFile << endl;

    // P(TXE)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] Total Excitation Energy Probability P(TXE)" << endl << endl;
    sum = 0;
    for (int i=0; i<=numberEnergies; i++) {
	sum += initialTXE[i];
    }
    for (int i=0; i<=numberEnergies; i++) {
	initialTXE[i] /= sum;
	outputFile << energyGrid[i] << " " << initialTXE[i] << endl;
    }
    outputFile << endl;

    // <TKE>(A)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] A   <TKE>(A) <Q>(A)" << endl << endl;
    double avTKE_ = 0.;
    outputFile << fixed << setprecision(4) ;
    for (int i=0; i<NUMA; i++) {
	if (numberEventsPerMass[i]!=0) {
	    TKEvsMass[i] /= (double) numberEventsPerMass[i];
	    avTKE_+=TKEvsMass[i]*massYields[i];
	    sTKEvsMass[i]=sqrt(sTKEvsMass[i]/numberEventsPerMass[i]-TKEvsMass[i]*TKEvsMass[i]);
	    QfA[i] /= numberEventsPerMass[i];
	    outputFile << setw(4) <<  i << " " << setw(6) << TKEvsMass[i] << " " << setw(7) << right << sTKEvsMass[i] << " " << setw(6) << QfA[i] << " " << numberEventsPerMass[i] << endl;
	}
    }
    outputFile << endl;
    cout << " Average TKE: " << avTKE_/2. << endl;

    writeYieldsInFile ("testYieldTable.dat");

    // <TXE>(A)
    outputFile << "# [gnuplot " << ++gnuplotIndex << "] <TXE>=f(A) " << endl << endl;
    for (int i=0; i<NUMA; i++) {
	if (numberEventsPerMass[i]!=0) {
	    TXEvsMass[i] /= (double) numberEventsPerMass[i];
	    outputFile << i << " " << TXEvsMass[i] << endl;
	}
    }
    outputFile << endl;

    outputFile.close();

/* Uncomment for writeout of Y(A,TKE)
    //-- Write out Y(A,TKE) in separate file
    outputFile.open("YATKE.CGMF.out");
    outputFile << "# Y(A,TKE)\n";
    double sumYields;
    gnuplotIndex = -1;
    for (int j=0; j<NUMTKE; j++) {
	sumYields = 0.0;
	for (int i=0; i<NUMA; i++) {sumYields += YieldsATKE[i][j];}
	if (sumYields!=0.0) {
	    outputFile << "# [gnuplot " << ++gnuplotIndex << "] ; TKE = " << j << " " << sumYields << "\n\n";
	    for (int i=0; i<NUMA; i++) {
		YieldsATKE[i][j] /= sumYields; // Renormalization
		outputFile << i << " " << YieldsATKE[i][j] << "\n";
	    }
	    outputFile << endl;
	}
    }
*/

    outputFile.close();
    // Printing ------------------------------------------------------------------------------
    return;
}
// -----------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************
 * computeNeutronSeparationEnergies
 *------------------------------------------------------------------------------
 * Computes the 1 and 2-neutron separation energies from the Audi mass table
 * which has exp. compiled by Audi and theory from FRDM. Uncomment out the
 * lines for comparison of theory-theory and exp.-exp. only.
 ******************************************************************************/
void FissionFragments::computeNeutronSeparationEnergies(void) {

    double Mn = mass_excess(0,1); // neutron mass excess (MeV)

    for (int i=Amin; i<=Amax; i++) {
	for (int j=Zmin; j<=Zmax; j++) {
	    if (YZA2[j][i]!=0.0) {
		S1n[j][i] = -mass_excess(j,i) + mass_excess(j,i - 1) + Mn;
		S2n[j][i] = -mass_excess(j,i) + mass_excess(j,i - 2) + 2*Mn;
	    }
	}
    }
    return;
}

/*******************************************************************************
 *                       SUBROUTINE computePfactors
 *------------------------------------------------------------------------------
 * Computes Pfactors for every fission fragment yield, and returns a mass-
 * dependent array Pf(Amin:Ac/2) that contains the Y(Z|A)-weighted Pfactors.
 ******************************************************************************/
void FissionFragments::computePfactors (void) {
    int Zpl, Zph;
    for (int i=Amin; i<=Asym; i++) {
	Zpl = (int) (Z0[i] + 0.5);
	Zph = Zc - Zpl;
	for (int j=-dZ; j<=dZ; j++) {
	    Pfactors[i][Zpl + j]      = Pfactor(Zpl + j, i - Zpl - j);
	    Pfactors[Ac - i][Zph - j] = Pfactor(Zph - j, Ac - i - Zph + j);
	}
    }

    ofstream out;
    out.open("Pfactors.out");
    double pfl, pfh, sum;
    double y;
    for (int i=Amin; i<=Asym; i++) {
	Zpl = (int) (Z0[i] + 0.5);
	Zph = Zc - Zpl;
	pfl = 0.0;
	pfh = 0.0;
	sum = 0.0;
	for (int j=-dZ; j<=dZ; j++) {
	    y = YZA2[Zpl + j][i];
	    pfl = pfl + y*Pfactors[i][Zpl + j];
	    pfh = pfh + y*Pfactors[Ac - i][Zph - j];
	    sum += y;
	}
	if (sum!=0.0) { pfl /= sum; pfh /= sum; }

    //    pfl = Pfactors[i][Zpl];
    //    pfh = Pfactors[Ac-i][Zc-Zpl];

    }
    out.close();
    return;
}

/*******************************************************************************
 * computeTemperatures
 *------------------------------------------------------------------------------
 * Computes temperatures for a given nucleus (Z,A) on a given energy grid
 * (with NUME number of energy points).
 * The temperature is calculated as: 1/T = ( d/dU ) log( LevelDensity( U )
 ******************************************************************************/
void FissionFragments::computeTemperatures (int Z, int A) {

//    cout << "# Computing Temperatures for " << Z << " " << A << "\n\n";

    Nucleus nucleus;

    nucleus.za.setZA(Z,A);
//    nucleus.ndisc = riplReadDiscreteLevels(&nucleus.za,nucleus.lev,reassign);
    getRiplDiscreteLevels (&nucleus);
    statFixDiscreteLevels(&nucleus);
    statSetupEnergyBin(&nucleus);
    statSetupLevelDensityParameter(&nucleus,&nucleus.ldp);

    double *ld;
    ld = new double [NUME];

    for (int k=0; k<NUME; k++) {
	ld[k] = ldLevelDensity(ldEnergyGrid[k], A, &nucleus.ldp);
    }

    double de;
    double *temp;
    temp = new double [NUME];

    double *ldp = new double [NUME];

    for (int k=0; k<NUME; k++) {
	ldp[k] = ldDensityParameter(ldEnergyGrid[k], A, &nucleus.ldp);
    }

    for (int k=2; k<NUME-2; k++) {
	de = ldEnergyGrid[k] - ldEnergyGrid[k-1];
	temp[k] = ( -log(ld[k+2]) + 8.0*log(ld[k+1]) - 8.0*log(ld[k-1]) + log(ld[k-2]) ) / (12.*de);
    }

    de = ldEnergyGrid[1]-ldEnergyGrid[0];
    temp[0] = ( -25.0*log(ld[0]) + 48.0*log(ld[1]) - 36.0*log(ld[2]) + 16.0*log(ld[3]) - 3.0*log(ld[4]) ) / (12.0*de);
  
    de = ldEnergyGrid[2]-ldEnergyGrid[1];
    temp[1] = ( -25.0*log(ld[1]) + 48.0*log(ld[2]) - 36.0*log(ld[3]) + 16.0*log(ld[4]) - 3.0*log(ld[5]) ) / (12.0*de);
  
    de = ldEnergyGrid[NUME-2]-ldEnergyGrid[NUME-3];
    temp[NUME-2] = ( 25.0*log(ld[NUME-2]) - 48.0*log(ld[NUME-3]) + 36.0*log(ld[NUME-4]) - 16.0*log(ld[NUME-5]) + 3.0*log(ld[NUME-6]) ) / (12.0*de);
  
    de = ldEnergyGrid[NUME-1]-ldEnergyGrid[NUME-2];
    temp[NUME-1] = ( 25.0*log(ld[NUME-1]) - 48.0*log(ld[NUME-2]) + 36.0*log(ld[NUME-3]) - 16.0*log(ld[NUME-4]) + 3.0*log(ld[NUME-5]) ) / (12.0*de);
  
    for (int k=0; k<NUME; k++) temp[k] = 1.0/temp[k];

/*  for (int k=0; k<NUME; k++) {
		cout << k << " " << ldEnergyGrid[k] << " " << ld[k] << " " << temp[k] << "\n";
  }
 */
    cout << 1000*Z+A << " ";

	/*	for (int k=0; k<NUME; k++) {
		cout << fixed << setprecision(2) << temp[k] << " ";
	}
	cout << endl;
*/
    for (int k=0; k<41; k+=4) { // 0-10 MeV per step of 1.0 MeV
	cout << fixed << setprecision(2) << ldp[k] << " ";
    }
    for (int k=48; k<81; k+=8) { // 11-20 MeV per step of 2.0 MeV
	cout << fixed << setprecision(2) << ldp[k] << " ";
    }
    for (int k=100; k<NUME; k+=20) { // 25-100 MeV per step of 5.0 MeV
	cout << fixed << setprecision(2) << ldp[k] << " ";
    }

    cout << endl;


    delete ld;
    delete temp;
    return;
}

/*******************************************************************************
 * computeTemperatureTable
 *------------------------------------------------------------------------------
 * Computes the temperature tables T(A,Z,U) for the fragments (Z,A) on a given
 * energy grid with NUME points. These are used in calculating the spin
 * parameter B through 1/T = (d/dU) log [ rho(A,Z,U) ]
 ******************************************************************************/
void FissionFragments::computeTemperatureTables(void) {
   
   cout << "[computeTemperatureTables] Running... ";
   
   double de;
   
   temperatures = new double **[NUMA];
   for (int i=Amin; i<=Amax; i++) {
      temperatures[i] = new double *[NUMdZ];
      for (int j=0; j<=2*dZ1; j++) {
         temperatures[i][j] = new double [NUME];
         //      if (YZA[j][i]!=0.0) {
         
         for (int k=2; k<NUME-2; k++) {
            de = ldEnergyGrid[k] - ldEnergyGrid[k-1];
            temperatures[i][j][k] =
            ( -log(levelDensities[i][j][k+2]) + 8.0*log(levelDensities[i][j][k+1]) -
             8.0*log(levelDensities[i][j][k-1]) + log(levelDensities[i][j][k-2]) ) / (12.*de);
            //          cout << " ... " << k << "  " << beta_temperature[i][j][k] << "\n";
         }
         
         de = ldEnergyGrid[1]-ldEnergyGrid[0];
         temperatures[i][j][0] =
         ( -25.0*log(levelDensities[i][j][0]) + 48.0*log(levelDensities[i][j][1]) -
          36.0*log(levelDensities[i][j][2]) + 16.0*log(levelDensities[i][j][3]) -
          3.0*log(levelDensities[i][j][4]) ) / (12.0*de);
         
         de = ldEnergyGrid[2]-ldEnergyGrid[1];
         temperatures[i][j][1] =
         ( -25.0*log(levelDensities[i][j][1]) + 48.0*log(levelDensities[i][j][2]) -
          36.0*log(levelDensities[i][j][3]) + 16.0*log(levelDensities[i][j][4]) -
          3.0*log(levelDensities[i][j][5]) ) / (12.0*de);
         
         de = ldEnergyGrid[NUME-2]-ldEnergyGrid[NUME-3];
         temperatures[i][j][NUME-2] =
         ( 25.0*log(levelDensities[i][j][NUME-2]) - 48.0*log(levelDensities[i][j][NUME-3]) +
          36.0*log(levelDensities[i][j][NUME-4]) - 16.0*log(levelDensities[i][j][NUME-5]) +
          3.0*log(levelDensities[i][j][NUME-6]) ) / (12.0*de);
         
         de = ldEnergyGrid[NUME-1]-ldEnergyGrid[NUME-2];
         temperatures[i][j][NUME-1] =
         ( 25.0*log(levelDensities[i][j][NUME-1]) - 48.0*log(levelDensities[i][j][NUME-2]) +
          36.0*log(levelDensities[i][j][NUME-3]) - 16.0*log(levelDensities[i][j][NUME-4]) +
          3.0*log(levelDensities[i][j][NUME-5]) ) / (12.0*de);
         
         for (int k=0; k<NUME; k++) {
            temperatures[i][j][k] = 1.0/temperatures[i][j][k];
            //	  cout << i << " " << j << " " << k << temperatures[i][j][k];
         }
         
         //      }
      }
   }
   
   cout << "OK" << endl;
   return;
}

/*******************************************************************************
 * Constructor for charge distribution
 *------------------------------------------------------------------------------
 * Initializes a charge distribution class, where the resulting vectors Z0, sZ0,
 * FZZ0, FNZ0 are set to 0. We initialize similar vectors internal to the class,
 * define the energy regions, and set the parameters from Wahl's systematics.
 ******************************************************************************/
ZDistribution::ZDistribution(double* Z0, double* sZ0, double* FZZ0, double* FNZ0) {

	// Reset the Zp parameters with zeroes
	std::fill_n(Z0,NUMA,0.0);
	std::fill_n(sZ0,NUMA,0.0);
	std::fill_n(FZZ0,NUMA,0.0);
	std::fill_n(FNZ0,NUMA,0.0);

	// Initialize and reset the temporary Zp parameters
	std::fill_n(delz,NUMA,0.0);
	std::fill_n(sigz,NUMA,0.0);
	std::fill_n(Fz,NUMA,0.0);
	std::fill_n(Fn,NUMA,0.0);

	// Energy boundaries between low, intermediate, and high-energy
	PE_L = 8.0; // 8 MeV
	PE_H = 20.0; // 20 MeV

	/* LOW-ENERGY REGION -- see Tab. 2 of Wahl's systematics */
	// Peak parameters
	sigz140_Lc.push_back( 0.566); sigz140_Lc.push_back(0.0); sigz140_Lc.push_back( 0.0064); sigz140_Lc.push_back(0.0109); sigz140_Lc.push_back( 0.0);
	delz140_Lc.push_back(-0.487); delz140_Lc.push_back(0.0); delz140_Lc.push_back( 0.0180); delz140_Lc.push_back(0.0);    delz140_Lc.push_back(-0.00203);

	// Legacy CGMF values
	  Fz140_Lc.push_back( 1.242);   Fz140_Lc.push_back(0.0);   Fz140_Lc.push_back(-0.0183);   Fz140_Lc.push_back(-0.0152);      Fz140_Lc.push_back( 0.0);
	// Adjust for the method - what is in Wahl's table
	if (Method == "NewCGMF") {
		Fz140_Lc[0] = 1.207; Fz140_Lc[1] = 0.0; Fz140_Lc[2] = -0.0420; Fz140_Lc[3] = 0.0; Fz140_Lc[4] = 0.0022;
	}

	  Fn140_Lc.push_back( 1.076);   Fn140_Lc.push_back(0.0);   Fn140_Lc.push_back( 0.0);      Fn140_Lc.push_back(0.0);      Fn140_Lc.push_back( 0.0);
	 sigzSL_Lc.push_back(-0.0038); sigzSL_Lc.push_back(0.0);  sigzSL_Lc.push_back( 0.0);     sigzSL_Lc.push_back(0.0);     sigzSL_Lc.push_back( 0.0);
	 delzSL_Lc.push_back(-0.0080); delzSL_Lc.push_back(0.0);  delzSL_Lc.push_back( 0.0);     delzSL_Lc.push_back(0.0);     delzSL_Lc.push_back( 0.0);
	   FzSL_Lc.push_back( 0.0030);   FzSL_Lc.push_back(0.0);    FzSL_Lc.push_back( 0.0);       FzSL_Lc.push_back(0.0);       FzSL_Lc.push_back( 0.0);
	   FnSL_Lc.push_back(-0.0006);   FnSL_Lc.push_back(0.0);    FnSL_Lc.push_back( 0.0);       FnSL_Lc.push_back(0.0);       FnSL_Lc.push_back( 0.0);

	// Symmetry parameters
	   SL50_Lc.push_back(0.191);    SL50_Lc.push_back(0.0);     SL50_Lc.push_back(-0.0076); SL50_Lc.push_back(0.0);    SL50_Lc.push_back(0.0);
	 sigz50_Lc.push_back(0.356);  sigz50_Lc.push_back(0.060); sigz50_Lc.push_back( 0.0);  sigz50_Lc.push_back(0.0);  sigz50_Lc.push_back(0.0);
	delzmax_Lc.push_back(0.699); delzmax_Lc.push_back(0.0);  delzmax_Lc.push_back( 0.0); delzmax_Lc.push_back(0.0); delzmax_Lc.push_back(0.0);

	// Wing parameters
	sigzSLW_Lc.push_back(-0.045); sigzSLW_Lc.push_back( 0.0094); sigzSLW_Lc.push_back(0.0); sigzSLW_Lc.push_back(0.0); sigzSLW_Lc.push_back(0.0);
	delzSLW_Lc.push_back( 0.0);   delzSLW_Lc.push_back(-0.0045); delzSLW_Lc.push_back(0.0); delzSLW_Lc.push_back(0.0); delzSLW_Lc.push_back(0.0);
	  FzSLW_Lc.push_back( 0.159);   FzSLW_Lc.push_back(-0.028);    FzSLW_Lc.push_back(0.0);   FzSLW_Lc.push_back(0.0);   FzSLW_Lc.push_back(0.0);
	  FnSLW_Lc.push_back( 0.039);   FnSLW_Lc.push_back( 0.0);      FnSLW_Lc.push_back(0.0);   FnSLW_Lc.push_back(0.0);   FnSLW_Lc.push_back(0.0);

	/* HIGH-ENERGY REGION -- see Tab. 3 of Wahl's systematics */
	// Peak parameters
	sigz140_Hc.push_back( 0.542); sigz140_Hc.push_back(1.310); sigz140_Hc.push_back(0.033); sigz140_Hc.push_back(0.0);   sigz140_Hc.push_back(-0.005);
	delz140_Hc.push_back(-0.428); delz140_Hc.push_back(0.0);   delz140_Hc.push_back(0.0);   delz140_Hc.push_back(0.164); delz140_Hc.push_back(-0.0116);
	  Fz140_Hc.push_back( 1.0);     Fz140_Hc.push_back(0.0);     Fz140_Hc.push_back(0.0);     Fz140_Hc.push_back(0.0);     Fz140_Hc.push_back( 0.0);
	  Fn140_Hc.push_back( 1.0);     Fn140_Hc.push_back(0.0);     Fn140_Hc.push_back(0.0);     Fn140_Hc.push_back(0.0);     Fn140_Hc.push_back( 0.0);
	 sigzSL_Hc.push_back( 0.0);    sigzSL_Hc.push_back(0.0);    sigzSL_Hc.push_back(0.0);    sigzSL_Hc.push_back(0.0);    sigzSL_Hc.push_back( 0.0);
	 delzSL_Hc.push_back( 0.0);    delzSL_Hc.push_back(0.0);    delzSL_Hc.push_back(0.0);    delzSL_Hc.push_back(0.0);    delzSL_Hc.push_back( 0.0);
	   FzSL_Hc.push_back( 0.0);      FzSL_Hc.push_back(0.0);      FzSL_Hc.push_back(0.0);      FzSL_Hc.push_back(0.0);      FzSL_Hc.push_back( 0.0);
	   FnSL_Hc.push_back( 0.0);      FnSL_Hc.push_back(0.0);      FnSL_Hc.push_back(0.0);      FnSL_Hc.push_back(0.0);      FnSL_Hc.push_back( 0.0);

	// Symmetry parameters
	   SL50_Hc.push_back(0.191);   SL50_Hc.push_back(0.0);     SL50_Hc.push_back(-0.0076);  SL50_Hc.push_back(0.0);    SL50_Hc.push_back( 0.0);
	 sigz50_Hc.push_back(0.542); sigz50_Hc.push_back(1.310); sigz50_Hc.push_back( 0.033); sigz50_Hc.push_back(0.0);  sigz50_Hc.push_back(-0.005);
	delzmax_Hc.push_back(0.0);  delzmax_Hc.push_back(0.0);  delzmax_Hc.push_back( 0.0);  delzmax_Hc.push_back(0.0); delzmax_Hc.push_back( 0.0);

	// Wing parameters
	sigzSLW_Hc.push_back(0.0); sigzSLW_Hc.push_back(0.0); sigzSLW_Hc.push_back(0.0); sigzSLW_Hc.push_back(0.0); sigzSLW_Hc.push_back(0.0);
	delzSLW_Hc.push_back(0.0); delzSLW_Hc.push_back(0.0); delzSLW_Hc.push_back(0.0); delzSLW_Hc.push_back(0.0); delzSLW_Hc.push_back(0.0);
	  FzSLW_Hc.push_back(0.0);   FzSLW_Hc.push_back(0.0);   FzSLW_Hc.push_back(0.0);   FzSLW_Hc.push_back(0.0);   FzSLW_Hc.push_back(0.0);
	  FnSLW_Hc.push_back(0.0);   FnSLW_Hc.push_back(0.0);   FnSLW_Hc.push_back(0.0);   FnSLW_Hc.push_back(0.0);   FnSLW_Hc.push_back(0.0);

	return;
}

/*******************************************************************************
 * Computes the parameters needed for a charge distribution
 *------------------------------------------------------------------------------
 * Given a charge Zf and mass Af of the fissioning nucleus, along with its
 * excitation energy PE, we determine what energy region this falls into. For
 * PE < PE_L, we use the low-energy parameters. For PE_L <= PE <= PE_H, we form
 * a linear interpolation between low and high-energy parameters. For PE > PE_H
 * we use the high-energy parameters.
 ******************************************************************************/
void ZDistribution::ComputeZpParams(int Zf, int Af, double PE) {

	// Low-energy fissioning
	if (PE <= PE_L) {
		sigz140 = WahlEq17(Zf,Af,PE,sigz140_Lc);
		delz140 = WahlEq17(Zf,Af,PE,delz140_Lc);
		Fz140   = WahlEq17(Zf,Af,PE,Fz140_Lc);
		Fn140   = WahlEq17(Zf,Af,PE,Fn140_Lc);
		sigzSL  = WahlEq17(Zf,Af,PE,sigzSL_Lc);
		delzSL  = WahlEq17(Zf,Af,PE,delzSL_Lc);
		FzSL    = WahlEq17(Zf,Af,PE,FzSL_Lc);
		FnSL    = WahlEq17(Zf,Af,PE,FnSL_Lc);
		SL50    = WahlEq17(Zf,Af,PE,SL50_Lc);
		sigz50  = WahlEq17(Zf,Af,PE,sigz50_Lc);
		delzmax = WahlEq17(Zf,Af,PE,delzmax_Lc);
		sigzSLW = WahlEq17(Zf,Af,PE,sigzSLW_Lc);
		delzSLW = WahlEq17(Zf,Af,PE,delzSLW_Lc);
		FzSLW   = WahlEq17(Zf,Af,PE,FzSLW_Lc);
		FnSLW   = WahlEq17(Zf,Af,PE,FnSLW_Lc);
	}
	// Intermediate energy fission
	if ((PE >= PE_L) && (PE <= PE_H)) {

		// Do a linear interpolation between low and high values
		double FRH = (PE - PE_L)/(PE_H - PE_L);
		double FRL = 1.0 - FRH;

		// The low-energy part
		sigz140_L = WahlEq17(Zf,Af,PE_L,sigz140_Lc);
		delz140_L = WahlEq17(Zf,Af,PE_L,delz140_Lc);
		Fz140_L   = WahlEq17(Zf,Af,PE_L,Fz140_Lc);
		Fn140_L   = WahlEq17(Zf,Af,PE_L,Fn140_Lc);
		sigzSL_L  = WahlEq17(Zf,Af,PE_L,sigzSL_Lc);
		delzSL_L  = WahlEq17(Zf,Af,PE_L,delzSL_Lc);
		FzSL_L    = WahlEq17(Zf,Af,PE_L,FzSL_Lc);
		FnSL_L    = WahlEq17(Zf,Af,PE_L,FnSL_Lc);
		SL50_L    = WahlEq17(Zf,Af,PE_L,SL50_Lc);
		sigz50_L  = WahlEq17(Zf,Af,PE_L,sigz50_Lc);
		delzmax_L = WahlEq17(Zf,Af,PE_L,delzmax_Lc);
		sigzSLW_L = WahlEq17(Zf,Af,PE_L,sigzSLW_Lc);
		delzSLW_L = WahlEq17(Zf,Af,PE_L,delzSLW_Lc);
		FzSLW_L   = WahlEq17(Zf,Af,PE_L,FzSLW_Lc);
		FnSLW_L   = WahlEq17(Zf,Af,PE_L,FnSLW_Lc);

		// The high-energy part
		sigz140_H = WahlEq19(Zf,Af,PE_H,sigz140_Hc);
		delz140_H = WahlEq19(Zf,Af,PE_H,delz140_Hc);
		Fz140_H   = WahlEq19(Zf,Af,PE_H,Fz140_Hc);
		Fn140_H   = WahlEq19(Zf,Af,PE_H,Fn140_Hc);
		sigzSL_H  = WahlEq19(Zf,Af,PE_H,sigzSL_Hc);
		delzSL_H  = WahlEq19(Zf,Af,PE_H,delzSL_Hc);
		FzSL_H    = WahlEq19(Zf,Af,PE_H,FzSL_Hc);
		FnSL_H    = WahlEq17(Zf,Af,PE_H,FnSL_Hc);
		SL50_H    = WahlEq19(Zf,Af,PE_H,SL50_Hc);
		sigz50_H  = WahlEq19(Zf,Af,PE_H,sigz50_Hc);
		delzmax_H = WahlEq19(Zf,Af,PE_H,delzmax_Hc);
		sigzSLW_H = WahlEq19(Zf,Af,PE_H,sigzSLW_Hc);
		delzSLW_H = WahlEq19(Zf,Af,PE_H,delzSLW_Hc);
		FzSLW_H   = WahlEq19(Zf,Af,PE_H,FzSLW_Hc);
		FnSLW_H   = WahlEq19(Zf,Af,PE_H,FnSLW_Hc);

		// Calculate the interpolated values
		sigz140 = FRL*sigz140_L + FRH*sigz140_H;
		delz140 = FRL*delz140_L + FRH*delz140_H;
		Fz140   = FRL*Fz140_L + FRH*Fz140_H;
		Fn140   = FRL*Fn140_L + FRH*Fn140_H;
		sigzSL  = FRL*sigzSL_L + FRH*sigzSL_H;
		delzSL  = FRL*delzSL_L + FRH*delzSL_H;
		FzSL    = FRL*FzSL_L + FRH*FzSL_H;
		FnSL    = FRL*FnSL_L + FRH*FnSL_H;
		SL50    = FRL*SL50_L + FRH*SL50_H;
		sigz50  = FRL*sigz50_L + FRH*sigz50_H;
		delzmax = FRL*delzmax_L + FRH*delzmax_H;
		sigzSLW = FRL*sigzSLW_L + FRH*sigzSLW_H;
		delzSLW = FRL*delzSLW_L + FRH*delzSLW_H;
		FzSLW   = FRL*FzSLW_L + FRH*FzSLW_H;
		FnSLW   = FRL*FnSLW_L + FRH*FnSLW_H;

		// In LegacyCGMF, there is no linear interpolation, just a sharp change to intermediate energy treatment
		if (Method == "LegacyCGMF") {
			sigz140 = WahlEq19(Zf,Af,PE,sigz140_Hc);
			delz140 = WahlEq19(Zf,Af,PE,delz140_Hc);
			Fz140   = WahlEq19(Zf,Af,PE,Fz140_Hc);
			Fn140   = WahlEq19(Zf,Af,PE,Fn140_Hc);
			sigzSL  = WahlEq19(Zf,Af,PE,sigzSL_Hc);
			delzSL  = WahlEq19(Zf,Af,PE,delzSL_Hc);
			FzSL    = WahlEq19(Zf,Af,PE,FzSL_Hc);
			FnSL    = WahlEq17(Zf,Af,PE,FnSL_Hc);
			SL50    = WahlEq19(Zf,Af,PE,SL50_Hc);
			sigz50  = WahlEq19(Zf,Af,PE,sigz50_Hc);
			delzmax = WahlEq19(Zf,Af,PE,delzmax_Hc);
			sigzSLW = WahlEq19(Zf,Af,PE,sigzSLW_Hc);
			delzSLW = WahlEq19(Zf,Af,PE,delzSLW_Hc);
			FzSLW   = WahlEq19(Zf,Af,PE,FzSLW_Hc);
			FnSLW   = WahlEq19(Zf,Af,PE,FnSLW_Hc);
		}

	}
	// High-energy fission
	if (PE >= PE_H) {
		sigz140 = WahlEq19(Zf,Af,PE,sigz140_Hc);
		delz140 = WahlEq19(Zf,Af,PE,delz140_Hc);
		Fz140   = WahlEq19(Zf,Af,PE,Fz140_Hc);
		Fn140   = WahlEq19(Zf,Af,PE,Fn140_Hc);
		sigzSL  = WahlEq19(Zf,Af,PE,sigzSL_Hc);
		delzSL  = WahlEq19(Zf,Af,PE,delzSL_Hc);
		FzSL    = WahlEq19(Zf,Af,PE,FzSL_Hc);
		FnSL    = WahlEq19(Zf,Af,PE,FnSL_Hc);
		SL50    = WahlEq19(Zf,Af,PE,SL50_Hc);
		sigz50  = WahlEq19(Zf,Af,PE,sigz50_Hc);
		delzmax = WahlEq19(Zf,Af,PE,delzmax_Hc);
		sigzSLW = WahlEq19(Zf,Af,PE,sigzSLW_Hc);
		delzSLW = WahlEq19(Zf,Af,PE,delzSLW_Hc);
		FzSLW   = WahlEq19(Zf,Af,PE,FzSLW_Hc);
		FnSLW   = WahlEq19(Zf,Af,PE,FnSLW_Hc);

		// In LegacyCGMF, the delz140 parameter is negative of the calculated value
		if (Method == "LegacyCGMF") {
			sigz140 = WahlEq19(Zf,Af,PE,sigz140_Hc);
			delz140 = -WahlEq19(Zf,Af,PE,delz140_Hc);
			Fz140   = WahlEq19(Zf,Af,PE,Fz140_Hc);
			Fn140   = WahlEq19(Zf,Af,PE,Fn140_Hc);
			sigzSL  = WahlEq19(Zf,Af,PE,sigzSL_Hc);
			delzSL  = WahlEq19(Zf,Af,PE,delzSL_Hc);
			FzSL    = WahlEq19(Zf,Af,PE,FzSL_Hc);
			FnSL    = WahlEq17(Zf,Af,PE,FnSL_Hc);
			SL50    = WahlEq19(Zf,Af,PE,SL50_Hc);
			sigz50  = WahlEq19(Zf,Af,PE,sigz50_Hc);
			delzmax = WahlEq19(Zf,Af,PE,delzmax_Hc);
			sigzSLW = WahlEq19(Zf,Af,PE,sigzSLW_Hc);
			delzSLW = WahlEq19(Zf,Af,PE,delzSLW_Hc);
			FzSLW   = WahlEq19(Zf,Af,PE,FzSLW_Hc);
			FnSLW   = WahlEq19(Zf,Af,PE,FnSLW_Hc);
		}

	}

	// printf("-----\n");
	// printf("Some Wahl parameters\n");
	// printf("sigz140t=%12.6e\tdelz140t=%12.6e\tFz140t=%12.6e\tFn140t=%12.6e\n",sigz140,delz140,Fz140,Fn140);
	// printf("sigzSLt=%12.6e\tdelzSLt=%12.6e\tFzSLt=%12.6e\tSL50t=%12.6e\tsigz50t=%12.6e\n",sigzSL,delzSL,FzSL,SL50,sigz50);
	// printf("delzmaxt=%12.6e\tsigzSLWt=%12.6e\tdelzSLWt=%12.6e\tFzSLWt=%12.6e\tFnSLWt=%12.6e\n",delzmax,sigzSLW,delzSLW,FzSLW,FnSLW);
	// printf("-----\n");

	return;
}

/*******************************************************************************
 * Equation 17 from Wahl's systematics
 *------------------------------------------------------------------------------
 * Simply calculates the parameter value given the phenomenological equation 17
 * in Wahl's systematics. Takes in the fissioning nucleus' charge Zf, mass Af,
 * and excitation energy PE. Also takes in the coefficients corresponding to
 * whatever parameter we are calculating.
 ******************************************************************************/
double ZDistribution::WahlEq17(int Zf, int Af, double PE, vector<double> coeffs) {

	// Parameter
	double par = 0.0;

	// First check that the coefficients are of sufficient size
	if (coeffs.size() != 5) {
		printf(" -- WARNING: Wahl Eq17 requires 5 coefficients only!");
		return par;
	}

	// Double versions of Zf, Af
	double Zfd = double(Zf);
	double Afd = double(Af);

	// Eq. 17
	par = coeffs[0] + coeffs[1]*(Zfd - 92.) + coeffs[2]*(Afd - 236.) + coeffs[3]*(PE - 6.551) + coeffs[4]*(Afd - 236.)*(Afd - 236.);
	return par;
}

/*******************************************************************************
 * Equation 19 from Wahl's systematics
 *------------------------------------------------------------------------------
 * Simply calculates the parameter value given the phenomenological equation 19
 * in Wahl's systematics. Takes in the fissioning nucleus' charge Zf, mass Af,
 * and excitation energy PE. Also takes in the coefficients corresponding to
 * whatever parameter we are calculating.
 ******************************************************************************/
double ZDistribution::WahlEq19(int Zf, int Af, double PE, vector<double> coeffs) {

	// Parameter
	double par = 0.0;

	// First check that the coefficients are of sufficient size
	if (coeffs.size() != 5) {
		printf(" -- WARNING: Wahl Eq19 requires 5 coefficients only!");
		return par;
	}

	// Double versions of Zf
	double Zfd = double(Zf);

	// Eq. 19
	double P1 = coeffs[0] + coeffs[2]*(Zfd - 92.);
	double P2 = coeffs[1] + coeffs[3]*(Zfd - 92.);
	par = P1 + (P2 - P1)*(1.0 - exp(-coeffs[4]*PE));
	return par;
}

/*******************************************************************************
 * Computes the mass boundaries in the charge distribution
 *------------------------------------------------------------------------------
 * The Wahl systematics have 9 mass regions, which are symmetrical across the
 * symmetric mass. We calculate these based on the fissioning nucleus' charge
 * Zf, mass Af, and excitation energy PE. The result: mass boundaries as floats
 * (B1,B2,B3,Ba,Bb,B4,B5,B6) and integers (iB1,iB2,iB3,iBa,iBb,iB4,iB5,iB6).
 ******************************************************************************/
void ZDistribution::ComputeMassBounds(int Zf, int Af, double PE) {

	double dAf = (double)Af;
	double dZf = (double)Zf;

	// See bottom of pg. 23
	double F1 = (250. - dAf)/14.;
	// Check limits on F1
	if (F1<0.0) {
		F1 = 0.0;
	}
	if (F1>1.0) {
		F1 = 1.0;
	}
	// Form the primary fragment mass corresponding to the max yield
	double F2 = 1.0 - F1;
	double AK1 = 50.0*(dAf/dZf) - delzmax/SL50;
	double AK2 = (50.0 - delzmax)*(dAf/dZf);
	double Amax = F1*AK1 + F2*AK2;
	AKH = Amax;

	// See Eq. 10a - 10h
	B1 = 70.0;
	B2 = 77.0 + 0.036*(dAf - 236.0);

	Ba = dAf - Amax;
	Bb = Amax;
	/* -- This is the formula seen in Wahl's paper, but in his code cyfp.f90 it is the above - the paper appears to be incorrect
	Ba = Amax;
	Bb = dAf - Amax;
	*/

	/* To match the legacy CGMF, we had flipped delz140, but the mass bounds requires the original delz140 */
	double delz140_tmp = delz140;
	if (Method == "LegacyCGMF" && PE >= PE_H) {
		delz140_tmp = -delz140;
	}

	B4 = (delzmax - delz140_tmp + Amax*(SL50) + 140.0*delzSL)/(SL50 + delzSL);
	B3 = dAf - B4;
	B5 = dAf - B2;
	B6 = dAf - B1;

	// Compute the same boundaries but as integers (in Legacy CGMF)
	iB1 = 70;
	iB2 = (int)floor(77+0.036*(Af-236)+0.5);
	iBa = int(Af-Amax+0.5);
	iBb = int(Amax+0.5);
	iB4 = (int)floor((delzmax - delz140_tmp + Amax*SL50 + 140.0*delzSL)/(SL50 + delzSL) + 0.5);
	iB3 = Af - iB4;
	iB5 = Af - iB2;
	iB6 = Af - iB1;

	// Check boundaries
	if ((B1<B2) && (B2<B3) && (B3<Ba) && (Ba<Bb) && (Bb<B4) && (B4<B5) && (B5<B6)) {
		dAf = dAf;
	} else {
		printf(" -- WARNING: Wahl mass bounds have issues: %.1f %.1f %.1f %.1f %.1f %.1f %.1f %.1f\n",B1,B2,B3,Ba,Bb,B4,B5,B6);
	}

	// For debugging
	// printf("-----\n");
	// printf("Some Wahl parameters\n");
	// printf("B1=%d\tB2=%d\tB3=%d\tBa=%d\tBb=%d\tB4=%d\tB5=%d\tB6=%d\n",iB1,iB2,iB3,iBa,iBb,iB4,iB5,iB6);
	// printf("F1=%12.6e\tF2=%12.6e\tAK1=%12.6e\tAK2=%12.6e\tApmax=%12.6e\n",F1,F2,AK1,AK2,Amax);
	// printf("-----\n");

	return;
}

/*******************************************************************************
 * Computes the Wahl Zp parameters (delz,sigz,Fz,Fn) in the peak region (Eq. 11)
 *------------------------------------------------------------------------------
 * Loops through the Amin and Amax values, computing the Wahl Zp parameters in
 * the peak region (B2-B3 and B4-B5). Calculation is done either for the Legacy
 * version or the New version.
 ******************************************************************************/
void ZDistribution::PeakRegion(int Amin, int Amax, int Af) {

	double Afdbl = (double)Af;

	// Loop through all masses - must do this to agree with Legacy version
	for (int A=Amin;A<=Amax;A++) {
		double Adbl = (double)A;
		double Acomp = Afdbl - Adbl;

		// Legacy CGMF
		if ((A>=iB2) && (A<=iB3)) {
			delz[A] = -delz140 + delzSL*(A - (Af - 140));
			sigz[A] = sigz140 - sigzSL*(A - (Af - 140));
			Fz[A] = Fz140 - FzSL*(A - (Af - 140));
			Fn[A] = Fn140 - FzSL*(A - (Af - 140)); // Legacy CGMF calculations used FzSL as FnSL
		} else if ((A>=iB4) && (A<=iB5)) {
			delz[A] = delz140 + delzSL*(A - 140.);
			sigz[A] = sigz140 + sigzSL*(A - 140.);
			Fz[A] = Fz140 + FzSL*(A - 140);
			Fn[A] = Fn140 + FzSL*(A - 140);
		}
		if (Method == "LegacyCGMF") continue;
		// Current version of CGMF
		if ((Adbl>=B2) && (Adbl<=B3)) {
			delz[A] = -(delz140 + delzSL*(Acomp - 140.));
			sigz[A] = sigz140 + sigzSL*(Acomp - 140.);
			Fz[A] = Fz140 + FzSL*(Acomp - 140.);
			Fn[A] = Fn140 + FnSL*(Acomp - 140.);
		} else if ((Adbl>=B4) && (Adbl<=B5)) {
			delz[A] = delz140 + delzSL*(Adbl - 140.);
			sigz[A] = sigz140 + sigzSL*(Adbl - 140.);
			Fz[A] = Fz140 + FzSL*(Adbl - 140.);
			Fn[A] = Fn140 + FnSL*(Adbl - 140.); // While not explicitly given, FnSL is provided in the footnote (b) of Tab. 2
		}
	}
	return;
}

/*******************************************************************************
 * Computes the Wahl Zp parameters (delz,sigz,Fz,Fn) in the symmetry region (Eq. 12)
 *------------------------------------------------------------------------------
 * Loops through the Amin and Amax values, computing the Wahl Zp parameters in
 * the symmetry region (B3-Ba and Ba-Bb and Bb-B4). Calculation is done either
 * for the Legacy version or the New version.
 ******************************************************************************/
void ZDistribution::SymmetryRegion(int Amin, int Amax, int Af) {

	double Afdbl = (double)Af;

	// Loop through all masses - must do this to agree with Legacy version
	for (int A=Amin;A<=Amax;A++) {
		double Adbl = (double)A;
		double Acomp = Afdbl - Adbl;

		// Legacy CGMF
		if ((A>iB3) && (A<=iBa)) {
			delz[A] = delz[iB3] - SL50*(A - iB3);
			sigz[A] = sigz50;
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		} else if ((A>iBa) && (A<iBb)) {
			delz[A] = delz[iBa] + (A-iBa)*(2.0*delz[iBa])/(iBa-iBb);
			sigz[A] = sigz140 - sigzSL*(140 - iBb);
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		} else if ((A>=iBb) && (A<iB4)) {
			delz[A] = delz[iB4] + SL50*(iB4 - A);
			sigz[A] = sigz50;
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		}
		if (Method == "LegacyCGMF") continue;
		// Current version of CGMF
		if ((Adbl>=B3) && (Adbl<=Ba)) {
			delz[A] = -delzmax - SL50*(Adbl - Ba);
			sigz[A] = sigz50;
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		} else if ((Adbl>=Ba) && (Adbl<=Bb)) {
			delz[A] = -delzmax + (Adbl-Ba)*2.0*delzmax/(Bb - Ba);
			sigz[A] = sigz140 - sigzSL*(140. - Bb);
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		} else if ((Adbl>=Bb) && (Adbl<=B4)) {
			delz[A] = delzmax - SL50*(Adbl - Bb);
			sigz[A] = sigz50;
			Fz[A] = 1.0;
			Fn[A] = 1.0;
		}
	}
	return;
}

/*******************************************************************************
 * Computes the Wahl Zp parameters (delz,sigz,Fz,Fn) in the wing region (Eq. 13)
 *------------------------------------------------------------------------------
 * Loops through the Amin and Amax values, computing the Wahl Zp parameters in
 * the wing region (B1-B2 and B5-B6). Calculation is done either for the Legacy
 * version or the New version.
 ******************************************************************************/
void ZDistribution::WingRegion(int Amin, int Amax, int Af) {

	double Afdbl = (double)Af;

	// Loop through all masses - must do this to agree with Legacy version
	for (int A=Amin;A<=Amax;A++) {
		double Adbl = (double)A;
		double Acomp = Afdbl - Adbl;

		// Legacy CGMF
		if ((A>iB1) && (A<iB2)) {
			delz[A] = delz[iB2] + delzSLW*(iB2 - A);
			sigz[A] = sigz[iB5] + sigzSLW*(iB2 - A);
			Fz[A] = Fz140 + FzSLW*(iB2 - A);
			Fn[A] = Fn140 + FnSLW*(iB2 - A);
		} else if ((A>iB5) && (A<iB6)) {
			delz[A] = delz[iB5] - delzSLW*(A - iB5);
			sigz[A] = sigz[iB5] + sigzSLW*(A - iB5);
			Fz[A] = Fz140 + FzSLW*(A - iB5);
			Fn[A] = Fn140 + FnSLW*(A - iB5);
		}
		if (Method == "LegacyCGMF") continue;
		// Current version of CGMF
		if ((Adbl>=B1) && (Adbl<=B2)) {
			delz[A] = -delzB5 + delzSLW*(Acomp - B5);
			sigz[A] = sigzB5 + sigzSLW*(B2 - Adbl);
			Fz[A] = Fz140 + FzSLW*(B2 - Adbl);
			Fn[A] = Fn140 + FnSLW*(B2 - Adbl);
		} else if ((Adbl>=B5) && (Adbl<=B6)) {
			delz[A] = delzB5 - delzSLW*(Adbl - B5);
			sigz[A] = sigzB5 + sigzSLW*(Adbl - B5);
			Fz[A] = Fz140 + FzSLW*(Adbl - B5);
			Fn[A] = Fn140 + FnSLW*(Adbl - B5);
		}
	}
	return;
}

/*******************************************************************************
 * Computes the Wahl Zp parameters (delz,sigz,Fz,Fn) in the far wing region (Eq. 14)
 *------------------------------------------------------------------------------
 * Loops through the Amin and Amax values, computing the Wahl Zp parameters in
 * the far wing region (<B1 and >B6). Calculation is done either for the Legacy
 * version or the New version.
 ******************************************************************************/
void ZDistribution::FarWingRegion(int Amin, int Amax, int Af) {

	double Afdbl = (double)Af;

	// Loop through all masses - must do this to agree with Legacy version
	for (int A=Amin;A<=Amax;A++) {
		double Adbl = (double)A;
		double Acomp = Afdbl - Adbl;

		// Legacy CGMF
		if (A<=iB1) {
			delz[A] = delz[iB2];
			sigz[A] = sigz[iB5];
			Fz[A] = Fz140;
			Fn[A] = Fn140;
		} else if (A>=iB6) {
			delz[A] = delz[iB5];
			sigz[A] = sigz[iB5];
			Fz[A] = Fz140;
			Fn[A] = Fn140;
		}
		if (Method == "LegacyCGMF") continue;
		// Current version of CGMF
		if (Adbl<=B1) {
			delz[A] = -delzB5;
			sigz[A] = sigzB5;
			Fz[A] = Fz140;
			Fn[A] = Fn140;
		} else if (Adbl>=B6) {
			delz[A] = delzB5;
			sigz[A] = sigzB5;
			Fz[A] = Fz140;
			Fn[A] = Fn140;
		}
	}
	return;
}

/*******************************************************************************
 * Main computation engine for the charge distribution
 *------------------------------------------------------------------------------
 * Takes in the fissioning nucleus' charge Zf, mass Af, excitation energy PE,
 * along with the min. and max. fragment masses (Amin, Amax) and the number of
 * Z values (dZ) from Zp. Calculates the vectors Z0, sZ0, FZZ0, FNZ0 needed for
 * a charge distribution and updates the min. and max. Z value (Zmin, Zmax).
 ******************************************************************************/
void ZDistribution::ComputeZDistribution(int Zf, int Af, double PE, int Amin, int Amax, int dZ, double* Z0, double* sZ0, double* FZZ0, double* FNZ0, int* Zmin, int* Zmax) {

	// cout << " PE=" << PE <<endl;

	// Compute the parameters given the CN charge, mass, and excitation energy
	ComputeZpParams(Zf,Af,PE);

	// Compute the boundaries in primary fragment mass
	ComputeMassBounds(Zf,Af,PE);

	// For debugging
	// printf("%8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",B1,B2,B3,Ba,Bb,B4,B5,B6);

	// Check one of the Zp parameters
	if (sigz50 <= 0.0) {
		sigz50 = sigz140;
	}

	// Calculate the constants in Wahl's original code
	sigzB5 = sigz140 + sigzSL*(B5 - 140.);
	delzB5 = delz140 + delzSL*(B5 - 140.);

	// Calculate the Zp parameters for each primary fragment mass
	PeakRegion(Amin,Amax,Af);
	SymmetryRegion(Amin,Amax,Af);
	WingRegion(Amin,Amax,Af);
	FarWingRegion(Amin,Amax,Af);
	/* NOTE -- requires an A-loop in each function to match Legacy CGMF */

	// Calculate the mean charge Z0, width sZ0, and Z/N even-odd effects FZZ0, FNZ0 for each mass
	double x = double(Zf)/double(Af);
	for (int a=Amin;a<=Amax;a++) {
		double adbl = (double)a;
		Z0[a]   = x*adbl+delz[a];
		sZ0[a]  = sigz[a];
		FZZ0[a] = Fz[a];
		FNZ0[a] = Fn[a];
		// For debugging
		// printf("%d %8.4f %8.4f %8.4f %8.4f\n",a,delz[a],sigz[a],Fz[a],Fn[a]);
	}

	// Update the Zmin and Zmax
	*Zmin = (int) floor(Z0[Amin]+0.5)-dZ;
	*Zmax = (int) floor(Z0[Amax]+0.5)+dZ;

	// For debugging
	// printf("-----\n");
  	// printf("Wahl Parameters:\n");
	// for (int i=Amin; i<=Amax; i++) {
	// 	printf("%d\t%12.6e\t%12.6e\t%12.6e\t%12.6e\n",i,Z0[i],sZ0[i],FZZ0[i],FNZ0[i]);
	// }
	// printf("-----\n");

	return;
}

/*******************************************************************************
 * Returns the CGMF method (either Legacy or NewCGMF)
 *------------------------------------------------------------------------------
 * Simply returns the string corresponding to the CGMF method.
 ******************************************************************************/
string ZDistribution::ReturnMethod() {
	return Method;
}

/*******************************************************************************
 * Destructor for the charge distribution class
 *------------------------------------------------------------------------------
 * Clears memory of the vectors used to calculate the Wahl parameters.
 ******************************************************************************/
ZDistribution::~ZDistribution() {

	sigz140_Lc.clear();	sigz140_Hc.clear();
	delz140_Lc.clear();	delz140_Hc.clear();
	  Fz140_Lc.clear();	  Fz140_Hc.clear();
	  Fn140_Lc.clear();	  Fn140_Hc.clear();
	 sigzSL_Lc.clear();	 sigzSL_Hc.clear();
	 delzSL_Lc.clear();	 delzSL_Hc.clear();
	   FzSL_Lc.clear();	   FzSL_Hc.clear();
	   SL50_Lc.clear();	   SL50_Hc.clear();
	 sigz50_Lc.clear();	 sigz50_Hc.clear();
	delzmax_Lc.clear();	delzmax_Hc.clear();
	sigzSLW_Lc.clear();	sigzSLW_Hc.clear();
	delzSLW_Lc.clear();	delzSLW_Hc.clear();
	  FzSLW_Lc.clear();	  FzSLW_Hc.clear();
	  FnSLW_Lc.clear();	  FnSLW_Hc.clear();
}

// TODO: PROBABLY SHOULD ALSO RECODE THE BELOW TO MAKE SURE IT FOLLOWS A AND A' NOTATION (P.J.)
// -----------------------------------------------------------------------------------------------------------------------------
/*******************************************************************************
 * computeWahlParametersForChargeDistribution
 *------------------------------------------------------------------------------
 * Calculate the Wahl Parameters {Fn(A), Fz(A), delz(A), sigz(A)} required for
 * computing the charge distribution Y(Z|A) for any given fissioning isotope.
 *------------------------------------------------------------------------------
 * Reference:
 * A.C.Wahl, "Systematics of Fission-Product Yields", Los Alamos Technical
 * Report, LA-13928 (2002).
 *------------------------------------------------------------------------------
 * From Zeke Blain, blaine2@rpi.edu, F95 original version (March 2012).
 ******************************************************************************/
void FissionFragments::computeWahlParametersForChargeDistribution (void) {
  
  double sigz140t, delz140t, Fz140t, Fn140t, sigzSLt, delzSLt, FzSLt, SL50t,
  sigz50t, delzmaxt, sigzSLWt, delzSLWt, FzSLWt, FnSLWt;
  
  double delz[NUMA], sigz[NUMA], Fz[NUMA], Fn[NUMA];
  
  double PE;
  
  double sigz140[5] = {0.566, 0.0, 0.0064, 0.0109, 0.0};
  double delz140[5] = {-0.487, 0.0, 0.0180, 0.0, -0.00203};
  //	double Fz140[5]   = {1.207, 0.0, -0.0420, 0.0, 0.0022};
  double Fz140[5]   = {1.242, 0.0, -0.0183, -0.0152, 0.0};
  double Fn140[5]   = {1.076, 0.0, 0.0, 0.0, 0.0};
  double sigzSL[5]  = {-0.0038, 0.0, 0.0, 0.0, 0.0};
  double delzSL[5]  = {-0.0080, 0.0, 0.0, 0.0, 0.0};
  double FzSL[5]    = {0.0030, 0.0, 0.0, 0.0, 0.0};
  double SL50[5]    = {0.191, 0.0, -0.0076, 0.0, 0.0};
  double sigz50[5]  = {0.356, 0.060, 0.0, 0.0, 0.0};
  double delzmax[5] = {0.699, 0.0, 0.0, 0.0, 0.0};
  double sigzSLW[5] = {-.045, 0.0094, 0.0, 0.0, 0.0};
  double delzSLW[5] = {0.0, -0.0045, 0.0, 0.0, 0.0};
  double FzSLW[5]   = {0.159, -0.028, 0.0, 0.0, 0.0};
  double FnSLW[5]   = {0.039, 0.0, 0.0, 0.0, 0.0};
  
  std::fill_n(Z0,NUMA,0.0);
  std::fill_n(sZ0,NUMA,0.0);
  std::fill_n(FZZ0,NUMA,0.0);
  std::fill_n(FNZ0,NUMA,0.0);
  
  //============================================================================
  // TEMPORARY... FOR ANYTHING BUT n+U-235 FISSION
  //============================================================================
  /*  if (Ac!=236) {
   
   for (int i=Amin; i<=Amax; i++) {
   delz[i]  = 0.5; // dZ
   if(i>Asym) delz[i] = -0.5; // dZ
   sigz[i] = sigmaZ; // width of distribution
   Fz[i] = 1.0; // FZ - odd-even factor
   Fn[i] = 1.0; // FN - odd-even factor
   }
   
   //-- save Wahl parameters
   double x = float(Zt)/Ac;
   for (int i=Amin; i<=Amax; i++) {
   Z0[i]   = (x*i+delz[i]);
   sZ0[i]  = sigz[i];
   FZZ0[i] = Fz[i];
   FNZ0[i] = Fn[i];
   //      cout << i << " " << Z0[i] << " " << sZ0[i] << "\n";
   }
   
   //-- Find true (Zmin, Zmax)
   Zmin = (int) floor(Z0[Amin]+0.5)-dZ;
   Zmax = (int) floor(Z0[Amax]+0.5)+dZ;
   
   return;
   
   }
   */
  
  //  PE = incidentEnergy + B1n[Zc][Ac];
  PE = incidentEnergy + SnCompound;
  //  cout << "Einc=" << incidentEnergy << " PE=" << PE <<endl;

  // Calculate input parameter for any given fissioning system
  // >> SHOULD 6.551 BE REPLACED BY neutronSeparationEnergy? I THINK SO...
  // WHAT ABOUT (92,236)... IS IT JUST VALID FOR U235?
  if (PE <= 8.0) {
    
    sigz140t = sigz140[0]+sigz140[1]*(Zt-92)+sigz140[2]*(Ac-236)+sigz140[3]*(PE-6.551)+sigz140[4]*(Ac-236)*(Ac-236);
    delz140t = delz140[0]+delz140[1]*(Zt-92)+delz140[2]*(Ac-236)+delz140[3]*(PE-6.551)+delz140[4]*(Ac-236)*(Ac-236);
    Fz140t   = Fz140[0]+Fz140[1]*(Zt-92)+Fz140[2]*(Ac-236)+Fz140[3]*(PE-6.551)+Fz140[4]*(Ac-236)*(Ac-236);
    Fn140t   = Fn140[0]+Fn140[1]*(Zt-92)+Fn140[2]*(Ac-236)+Fn140[3]*(PE-6.551)+Fn140[4]*(Ac-236)*(Ac-236);
    sigzSLt  = sigzSL[0]+sigzSL[1]*(Zt-92)+sigzSL[2]*(Ac-236)+sigzSL[3]*(PE-6.551)+sigzSL[4]*(Ac-236)*(Ac-236);
    delzSLt  = delzSL[0]+delzSL[1]*(Zt-92)+delzSL[2]*(Ac-236)+delzSL[3]*(PE-6.551)+delzSL[4]*(Ac-236)*(Ac-236);
    FzSLt    = FzSL[0]+FzSL[1]*(Zt-92)+FzSL[2]*(Ac-236)+FzSL[3]*(PE-6.551)+FzSL[4]*(Ac-236)*(Ac-236);
    SL50t    = SL50[0]+SL50[1]*(Zt-92)+SL50[2]*(Ac-236)+SL50[3]*(PE-6.551)+SL50[4]*(Ac-236)*(Ac-236);
    sigz50t  = sigz50[0]+sigz50[1]*(Zt-92)+sigz50[2]*(Ac-236)+sigz50[3]*(PE-6.551)+sigz50[4]*(Ac-236)*(Ac-236);
    delzmaxt = delzmax[0]+delzmax[1]*(Zt-92)+delzmax[2]*(Ac-236)+delzmax[3]*(PE-6.551)+delzmax[4]*(Ac-236)*(Ac-236);
    sigzSLWt = sigzSLW[0]+sigzSLW[1]*(Zt-92)+sigzSLW[2]*(Ac-236)+sigzSLW[3]*(PE-6.551)+sigzSLW[4]*(Ac-236)*(Ac-236);
    delzSLWt = delzSLW[0]+delzSLW[1]*(Zt-92)+delzSLW[2]*(Ac-236)+delzSLW[3]*(PE-6.551)+delzSLW[4]*(Ac-236)*(Ac-236);
    FzSLWt   = FzSLW[0]+FzSLW[1]*(Zt-92)+FzSLW[2]*(Ac-236)+FzSLW[3]*(PE-6.551)+FzSLW[4]*(Ac-236)*(Ac-236);
    FnSLWt   = FnSLW[0]+FnSLW[1]*(Zt-92)+FnSLW[2]*(Ac-236)+FnSLW[3]*(PE-6.551)+FnSLW[4]*(Ac-236)*(Ac-236);
    
  } else if (PE<=30.) { // I Stetcu: this was extended to 30 MeV (see pages 22 and 25 for PE>~20 MeV)
    
    sigz140[0] = 0.542; sigz140[1] = 1.310; sigz140[2] = 0.033; sigz140[3] = 0.0; sigz140[4] = -0.005;
    delz140[0] = -0.428; delz140[1] = 0.0; delz140[2] = 0.0; delz140[3] = 0.164; delz140[4] = -0.0116;
    SL50[0] = 0.191; SL50[1] = 0.0; SL50[2] = -0.0076; SL50[3] = 0.0; SL50[4] = 0.0;
    sigz50[0] = 0.542; sigz50[1] = 1.310; sigz50[2] = 0.033; sigz50[3] = 0.0; sigz50[4] = -0.005;
    
    sigz140t = (sigz140[0]+sigz140[2]*(Zt-92))+((sigz140[1]+sigz140[3]*(Zt-92))-(sigz140[0]+sigz140[2]*(Zt-92)))*(1.0-exp(-sigz140[4]*PE));
    delz140t = (delz140[0]+delz140[2]*(Zt-92))+((delz140[1]+delz140[3]*(Zt-92))-(delz140[0]+delz140[2]*(Zt-92)))*(1.0-exp(-delz140[4]*PE));
    Fz140t   = 1.0;
    Fn140t   = 1.0;
    sigzSLt  = 0.0;
    delzSLt  = 0.0;
    FzSLt    = 0.0;
    SL50t    = (SL50[0]+SL50[2]*(Zt-92))+((SL50[1]+SL50[3]*(Zt-92))-(SL50[0]+SL50[2]*(Zt-92)))*(1.0-exp(-SL50[4]*PE));
    sigz50t  = (sigz50[0]+sigz50[2]*(Zt-92))+((sigz50[1]+sigz50[3]*(Zt-92))-(sigz50[0]+sigz50[2]*(Zt-92)))*(1.0-exp(-sigz50[4]*PE));
    delzmaxt = 0.0;
    sigzSLWt = 0.0;
    delzSLWt = 0.0;
    FzSLWt   = 0.0;
    FnSLWt   = 0.0;
    
    
  } else {

    cerr << "[computeWahlParameters] Wahl parameters only defined up to Einc+Sn=20 MeV!" << endl;
    exit(-1);
    
  }
  
  // Calculate the region values
  
  double F1=floor((250.-Ac)/14.+0.5);  // F1=anint((250.-Ac)/14.); in F95
  double F2=1.-F1;
  double AK1=50.0*float(Ac)/Zt-delzmaxt/SL50t;
  double AK2=(50.0-delzmaxt)*float(Ac)/Zt;
  double Apmax=F1*AK1+F2*AK2;
  
  int B1=70;
  int B2= (int) floor(77+0.036*(Ac-236)+0.5);      // nint() in F95
  int B4= (int) floor((delzmaxt-delz140t+Apmax*SL50t+140*delzSLt)/(SL50t+delzSLt)+0.5); // nint() in F95
  int B3=(Ac-B4);
  int B5=(Ac-B2);
  int B6=(Ac-B1);
  int Bb=int(Apmax+0.5); // nint() in F95
  int Ba=int(Ac-Apmax+0.5); // nint() in F95

  // Calculate values of sigz, delz, Fn, Fz for give Ap values
  for (int i=Amin; i<=Amax; i++) {
    // Peak Regions
    if ((i>=B2 && i<=B3) || (i>=B4 && i<=B5)) {
      if(PE<=20.){
	if (i>Ac/2.0) {
	  delz[i]=delz140t+delzSLt*(i-140.0);
	  sigz[i]=sigz140t+sigzSLt*(i-140.0);
	  Fz[i]=Fz140t+FzSLt*(i-140.0);
	  Fn[i]=Fn140t+FzSLt*(i-140.0);
	} else if (i<Ac/2.0) {
	  delz[i]=-1*delz140t+delzSLt*(i-(Ac-140));
	  sigz[i]=sigz140t-sigzSLt*(i-(Ac-140));
	  Fz[i]=Fz140t-FzSLt*(i-(Ac-140));
	  Fn[i]=Fn140t-FzSLt*(i-(Ac-140));
	}
      }else{
	sigz[i]=sigz140t;
	Fz[i]=1.;
	Fn[i]=1.;
	if(i<Ac/2.)
	  delz[i]=delz140t;
	else
	  delz[i]=-delz140t;
      }
    }
  }


  // Near the symmetry region, PE>20 uses the same parameterization as 8<PE<=20.  
  for (int i=Amin; i<=Amax; i++) {
    // Near Symmetry Region
    if (i>B3 && i<B4) {
      Fn[i]=1.;
      Fz[i]=1.;
      if (i>B3 && i<=Ba) {
        delz[i]=delz[B3]-SL50t*(i-B3);
        sigz[i]=sigz50t;
      } else if (i>Ba && i<Bb) {
        delz[i]=delz[Ba]+(i-Ba)*(2.0*delz[Ba])/(Ba-Bb);
        sigz[i]=sigz140t-sigzSLt*(140-Bb);
      } else if (i>=Bb && i<B4) {
        delz[i]=delz[B4]+SL50t*(B4-i);
        sigz[i]=sigz50t;
      }
    }
  }
  
  for (int i=Amin; i<=Amax; i++) {
    // Wing Regions

    if ((i>B1 && i<B2) || (i>B5 && i<B6)) {
      if(PE<=20.){
	if (i>Ac/2.0) {
	  delz[i]=delz[B5]-delzSLWt*(i-B5);
	  sigz[i]=sigz[B5]+sigzSLWt*(i-B5);
	  Fz[i]=Fz140t+FzSLWt*(i-B5);
	  Fn[i]=Fn140t+FnSLWt*(i-B5);
	} else if (i<Ac/2.0) {
	  delz[i]=delz[B2]+delzSLWt*(B2-i);
	  sigz[i]=sigz[B5]+sigzSLWt*(B2-i);
	  Fz[i]=Fz140t+FzSLWt*(B2-i);
	  Fn[i]=Fn140t+FnSLWt*(B2-i);
	}
      }else{
	sigz[i]=sigz140t;
	Fz[i]=1.;
	Fn[i]=1.;
	if(i<Ac/2.)
	  delz[i]=delz140t;
	else
	  delz[i]=-delz140t;
      }
      
    }

    
    // Far Wing Regions
    if (i<=B1 || i>=B6) {
      if(PE<=20.){
	if (i>Ac/2.0) {
	  delz[i]=delz[B5];
	  sigz[i]=sigz[B5];
	  Fz[i]=Fz140t;
	  Fn[i]=Fn140t;
	} else if (i<Ac/2.0) {
	  delz[i]=delz[B2];
	  sigz[i]=sigz[B5];
	  Fz[i]=Fz140t;
	  Fn[i]=Fn140t;
	}
      }else{
	sigz[i]=sigz140t;
	Fz[i]=1.;
	Fn[i]=1.;
	if(i<Ac/2.)
	  delz[i]=delz140t;
	else
	  delz[i]=-delz140t;
      }
    }
  }
  
  //-- save Wahl parameters
  double x = float(Zt)/Ac;
  for (int i=Amin; i<=Amax; i++) {
    Z0[i]   = (x*i+delz[i]);
    sZ0[i]  = sigz[i];
    FZZ0[i] = Fz[i];
    FNZ0[i] = Fn[i];
    //    FZZ0[i] = 1.0;
    //    FNZ0[i] = 1.0;
    //    cout << i << " " << Z0[i] << " " << sZ0[i] << "\n";
  }
  
  //-- Find true (Zmin, Zmax)
  Zmin = (int) floor(Z0[Amin]+0.5)-dZ;
  Zmax = (int) floor(Z0[Amax]+0.5)+dZ;
  
}
// -----------------------------------------------------------------------------------------------------------------------------

/*******************************************************************************
 * Pfactor
 *------------------------------------------------------------------------------
 * Computes the Pfactor for a nucleus of charge Z and neutron number N
 ******************************************************************************/
double FissionFragments::Pfactor(int Z, int N) {

    int Nn, Np; // Number of valence neutrons and protons
    int i;
    double p;

    const int magicNumbers[3] = {26,50,82}; // Spherical closed shells

    Nn = 99; Np = 99;
    // Solve for the number of valence nucleons
    for (int i=0; i<3; i++) {
	Nn = min(Nn,abs(N - magicNumbers[i]));
	Np = min(Np,abs(Z - magicNumbers[i]));
    }

    // Closed shell
    if (Np + Nn==0) {
	p = 0.0;
    } else { // Else it is Np*Nn/(Np+Nn)
	p = Np*Nn/float(Np+Nn);
    }

    return (p);
}


/*******************************************************************************
 * readDeformations
 *------------------------------------------------------------------------------
 * Read ground-state deformations (beta2 parameter) from FRDM95 calculations.
 * << * >> ADD REFERENCE!
 ******************************************************************************/
void readDeformations (void) {
	
    ifstream fp;
    string str = DEFORMATIONSFILE;
    int Z, A;

    // Try current directory
    fp.open(&str[0]);

    // Then system data area
    if (!fp) {
	str = datadir + str;
	fp.open(&str[0]);
    }

    if (!fp) cgmTerminateCode("deformtions data file not found");

    std::fill_n (beta2, MAX_ZAID, -100.0);

    while(getline(fp,str)){
	if(str[0] == '#') continue;
	Z = atoi(str.substr(0,4).c_str()); if (Z>99) break;
	A = atoi(str.substr(4,8).c_str()); if (A>199) continue;
	str = str.substr(44);
	if (str[6]!=' ') beta2[1000*Z+A] = atof(str.c_str());
    }

    fp.close();
    return;
}

/*******************************************************************************
 * getDeformation
 *------------------------------------------------------------------------------
 * Retrieve the beta2 parameter for a specific nucleus, identified by ZAID.
 ******************************************************************************/
double getDeformation (int zaid) {
	
    if (beta2[zaid]==-100.0) {
	cout << "Deformation parameter missing: " << zaid << endl;
	beta2[zaid]=0.2; // Default if missing
//	cgmTerminateCode ("No deformation parameter found for this nucleus: ");
    }
    return beta2[zaid];
}

/*******************************************************************************
 * readTemperatures
 *------------------------------------------------------------------------------
 * Read nuclear temperatures from a data file.
 * << * >> PROVIDE EXPLANATION ON THE ORIGIN OF THESE TEMPERATURES.
 ******************************************************************************/
void readTemperatures (void) {

    ifstream fp;
    string str = TEMPERATURESFILE;
    string s;
    int zaid;

    // Try current directory
    fp.open(&str[0]);

    // Then system data area
    if (!fp) {
	   str = datadir + str;
	   fp.open(&str[0]);
    }

    if (!fp) cgmTerminateCode("temperatures data file not found");

    for (int i=0; i<32; i++) {
  		for (int j=0; j<MAX_ZAID; j++) {
		    temperatures[j][i]=0.0;
		  }
    }

    int nc=0;
    while(getline(fp,str)){
      if(str[0] == '#') { nc++; } else { break; }
    }

    fp.clear();
    fp.seekg(0, ios::beg);

    for (int i=0; i<nc; i++) getline(fp,str);
    for (int i=0; i<NUM_GRID; i++) fp >> temperatureEnergyGrid[i];

    getline(fp,str); // skip comment line

    int i;
    while (getline(fp,str)) {
    	i=0;
    	while (!fp.eof()) {
  			fp >> zaid;
	   		for (int i=0; i<32; i++) fp >> temperatures[zaid][i];
    	}
    }

    return;
}

/*******************************************************************************
 * getTemperature
 *------------------------------------------------------------------------------
 * Retrieve the temperature for a specific nucleus at a given excitation energy.
 ******************************************************************************/
double getTemperature (int zaid, double energy) {

    if (temperatures[zaid][0]==0.0) {
	cout << "Temperature missing: " << zaid << endl;
	for (int i=0; i<NUM_GRID; i++) temperatures[zaid][i]=5.0; // Default if missing
//	cgmTerminateCode("temperature not found for this nucleus: ");
    }

    int iU = findEnergyIndex(energy, temperatureEnergyGrid);

    double e1 = temperatureEnergyGrid[iU];
    double e2 = temperatureEnergyGrid[iU+1];
    double t1 = temperatures[zaid][iU];
    double t2 = temperatures[zaid][iU+1];
    double t = (t2-t1)/(e2-e1)*(energy-e1)+t1;

    return t;
}

/*******************************************************************************
 * readSpinScaling
 *------------------------------------------------------------------------------
 * Reads the constant and linear coefficients for the alpha(En) dependence.
 * Values for alpha_0 and alpha_slope determined from fits to EgT(En) data for
 * U235(n,f) and Np237(n,f) data (Frehaut), or fits to Pnug.
 ******************************************************************************/
void readSpinScaling(void) {

	ifstream fp;
	string str = SPINSCALINGFILE;
  string line;
  int NumReactions = 0;

	// Try current directory
	fp.open(&str[0]);

	// The system data area
	if (!fp) {
		str = datadir + str;
		fp.open(&str[0]);
	}
	if (!fp) cgmTerminateCode("spin-scaling model file not found");

  while (getline(fp,line)) {
    if (line[0]=='#' or line=="") continue;
    istringstream(line) >> alpha_coefficients[NumReactions][0] >> alpha_coefficients[NumReactions][1] >> alpha_coefficients[NumReactions][2];
    NumReactions++;
    if (NumReactions>ALPHA_MAX_ZAIDc) cgmTerminateCode("Need to increase the size of alpha_coefficients array\n");
  }

	return;
}

/*******************************************************************************
 * readLDP
 *------------------------------------------------------------------------------
 * Read level density parameters from a data file.
 * << * >> PROVIDE EXPLANATION ON THE ORIGIN OF THESE TEMPERATURES.
 ******************************************************************************/
void readLDP (void) {

    ifstream fp;
    string str = LEVELDENSITYFILE;
    int zaid;

    // Try current directory
    fp.open(&str[0]);

    // Then system data area
    if (!fp) {
	str = datadir + str;
	fp.open(&str[0]);
    }

    if (!fp) cgmTerminateCode("level density parameters data file not found");

    for (int i=0; i<NUM_GRID; i++) {
		for (int j=0; j<MAX_ZAID; j++) {
		    ldp2[j][i]=0.0;
		}
    }

    while(getline(fp,str)){
		if(str[0] == '#') continue;
	    while (!fp.eof()) {
			fp >> zaid;
			for (int i=0; i<NUM_GRID; i++) fp >> ldp2[zaid][i];
    	}
	}

    return;
}

/*******************************************************************************
 * getLDP
 *------------------------------------------------------------------------------
 * Retrieve the level density parameter for a specific nucleus at a given
 * excitation energy.
 ******************************************************************************/
double getLDP (int zaid, double energy) {
	
    if (ldp2[zaid][0]==0.0) {
	cout << "ldp missing: " << zaid << endl;
	for (int i=0; i<NUM_GRID; i++) ldp2[zaid][i]=5.0; // Default if missing
//	cgmTerminateCode("level density parameter not found for this nucleus: ", zaid);
    }

    int iU = findEnergyIndex(energy, temperatureEnergyGrid);

    double e1 = temperatureEnergyGrid[iU];
    double e2 = temperatureEnergyGrid[iU+1];

    double p1 = ldp2[zaid][iU];
    double p2 = ldp2[zaid][iU+1];
    double p = (p2-p1)/(e2-e1)*(energy-e1)+p1;

    return p;
}

/*******************************************************************************
 * readNeutronBindingEnergies
 *------------------------------------------------------------------------------
 * Reads the 1-n and 2-n separation energies from the Audi-Wapstra tables. See
 * Ref. G. Audi et al., NPA 729 (2003) 337-676. Table contains exp. values.
 ******************************************************************************/
void FissionFragments::readNeutronBindingEnergies(void) {

    ifstream dataFile;
    string str;
    string string1, string2;

    string dataFilename = datadir;
    dataFilename.append("neutronSeparationEnergies.Audi2003.dat");

    int A, Z;
    char Symbol[2];
    double B1, B2;

    // Initialization
    for (int i=0; i<NUMA; i++) {
	for (int j=0; j<NUMZ; j++) {
	    B1n[j][i] = -99.0;
	    B2n[j][i] = -99.0;
	}
    }

    dataFile.open(&dataFilename[0]);
    // Check for file existence
    if (!dataFile) {cerr << "ERROR: Cannot find neutron separation energies from Audi 2003 file!" << endl; exit(-1);}

    while (!dataFile.eof()) { // Loop through file

	getline(dataFile, str);
	if (str.substr(0,1)!="#" && !str.empty()) { // Ignore comments or empty lines

	    istringstream(str) >> A >> Symbol >> Z; // Load A and Z

	    string1 = str.substr(12,9);
	    string2 = str.substr(32,9);

	    if (string1.find("*")==string::npos) { // Table contains "*" when there is no exp. value
		istringstream(string1) >> B1;
		B1n[Z][A] = B1 * 1e-3; // keV to MeV
	    }
	    if (string2.find("*")==string::npos) {
		istringstream(string2) >> B2;
		B2n[Z][A] = B2 * 1e-3; // keV to MeV
	    }
	}
    }
    dataFile.close();
    return;   
}

/*******************************************************************************
 * readAnisotropy
 *------------------------------------------------------------------------------
 * Read in the fission fragment anisotropy coefficients from anisotropy.dat
 ******************************************************************************/
void readAnisotropy (void) {

   ifstream fp;
   string line;
   int iZID = 1;
   int counter = 0;
   double Emin,Emax,Estep;
   int nEsteps;
   
   string str = ANISOTROPYFILE;
   // try to open in current directory
   fp.open(&str[0]);
   // otherwise, from data directory
   if (!fp) {
      str = datadir + str;
      fp.open(&str[0]);
   }


   while(getline(fp,line)){
      if(line[0] != '#') break;
   }

   istringstream(line) >> Emin >> Emax >> Estep;
   nEsteps = (int) (Emax-Emin)/Estep + 1; // compute number of energy steps

   for (int i=0;i<nEsteps;i++){
      anisotropy[0][i+1] = Estep*i;
   }

   while (!fp.eof()) {
      if (counter>=nEsteps) { 
	 iZID += 1;
	 counter = 0;
      }
      fp >> anisotropy[iZID][0];
      for (counter=0;counter<nEsteps;counter++){
         fp >> anisotropy[iZID][counter+1];
      }
   }

   return;

}

/*******************************************************************************
 * readPreEqAngularDistributionParameters
 *------------------------------------------------------------------------------
 * Read in the parameters for the fits to inelastic angular distributions
 ******************************************************************************/
void readPreEqAngularDistributionParameters (void) {

   ifstream fp;
   string line;

   string str = PREEQSCATTERFILE;
   // try to open in current directory
   fp.open(&str[0]);
   // otherwise, from data directory
   if (!fp) {
      str = datadir + str;
      fp.open(&str[0]);
   }

   while(getline(fp,str)){
   	if(str[0] == '#') continue;
   }

   // read in the parameters and store them into an array
   for(int i=0; i<4;i++){
      fp >> preeqScatteringParams[0][i][0] >> preeqScatteringParams[0][i][1] >> preeqScatteringParams[0][i][2] >> preeqScatteringParams[0][i][3] >> preeqScatteringParams[0][i][4] >> preeqScatteringParams[0][i][5];
   }
   for(int i=0; i<4;i++){
      fp >> preeqScatteringParams[1][i][0] >> preeqScatteringParams[1][i][1] >> preeqScatteringParams[1][i][2] >> preeqScatteringParams[1][i][3] >> preeqScatteringParams[1][i][4] >> preeqScatteringParams[1][i][5];
   }

   fp.close();

   return;

}

/*******************************************************************************
 * computeFragmentAngularDistribution
 *------------------------------------------------------------------------------
 * Selects the value of the anisotropy coefficient based on the energy of the 
 * incident neutron and the fissioning isotope
 * Calculates the P(cos(theta)) CDF to sample fragment angles
 * defined as P(cos(theta)) = norm*(1-ansiotropyCoefficient)*cos(theta)^2
 * norm is defined such that the sum over P(cos(theta))=1
 ******************************************************************************/
void FissionFragments::computeFragmentAngularDistribution (void) {

   // default value is 1.0 (isotropic)
   anisotropyCoefficient = 1.0;

   // determine the value for the fission fragment anisotropy (if tabulated)
   if (anisotropy_flag){
      for (int i=0;i<MAX_ZAIDc;i++){
         if (anisotropy[i+1][0]==ZAIDc-1){
            for (int j=1;j<=NUM_EGRID;j++){
               if (roundf(incidentEnergy0*10)/10 == anisotropy[0][j]) {
                  anisotropyCoefficient = anisotropy[i+1][j];
                  break;
               }
            }
         }
      }
   }

   // calculate the probability density function
   // if the anisotropy coefficient is not equal to 1.0
   if (anisotropyCoefficient != 1.0){
      double pdf = 0.0;
      double cosTheta;
      Pcostheta[0][0] = -1.0;
      Pcostheta[0][1] = 0.0;
      for (int i=1;i<THETA_STEPS;i++){
         cosTheta = 2.0*i/(THETA_STEPS-1)-1.0;
         pdf = 1.0+(anisotropyCoefficient-1.0)*cosTheta*cosTheta;
         Pcostheta[i][0] = cosTheta;
         Pcostheta[i][1] = Pcostheta[i-1][1]+pdf;
      }

      // normalize CDF
      for (int i=0;i<THETA_STEPS;i++){
         Pcostheta[i][1] = Pcostheta[i][1]/Pcostheta[THETA_STEPS-1][1];
      }
   }
   
   return ;

}

/*******************************************************************************
 * readMasses
 *------------------------------------------------------------------------------
 * Read excess masses from the Audi-Wapstra 2003 atomic mass evaluation.
 * Ref. G.Audi et al., NPA 729 (2003) 3-128.
 * << * >> IS IT AUDI 2011 NOW INSTEAD?
 * << * >> ADD REFERENCE!
 ******************************************************************************/
void readMasses (void) {

    ifstream fp;
    string d;
    int zaid;
    double m;

    string str = MASSESFILE;

    // Try current directory
    fp.open(&str[0]);

    // Then system data area
    if (!fp) {
	str = datadir + str;
	fp.open(&str[0]);
    }

    if (!fp) cgmTerminateCode("masses data file not found");

    std::fill_n (masses, MAX_ZAID, 0.0); // initialize to zeros

    while(getline(fp,str)){
	if(str[0] == '#') continue;
	while (str!="") {
	    zaid = atoi(str.substr(0,7).c_str());
	    m    = atof(str.substr(7,18).c_str());
	    masses[zaid] = m;
	    str = str.substr(18);
	}
    }

    return;
}

/*******************************************************************************
 * getMassExcess
 *------------------------------------------------------------------------------
 * Retrieve the mass excess (in MeV) for a specific nucleus.
 ******************************************************************************/
double getMassExcess (int zaid) {
  if (masses[zaid]==0.0){
    stringstream zaid_s;
    zaid_s<< zaid;
    cgmTerminateCode("Mass not found in mass table for this nucleus: "+zaid_s.str());
  }
    return masses[zaid];
}

/*******************************************************************************
 * readYieldsATKE
 *------------------------------------------------------------------------------
 * Reads the Y(A,TKE) from a file yieldsFilename
 ******************************************************************************/
void FissionFragments::readYieldsATKE(string yieldsFilename) {

    ifstream yieldsFile;
    int A, TKE;
    double yields;

    yieldsFile.open(&yieldsFilename[0]); // Check for file existence
    if (!yieldsFile) {cerr << "ERROR: Cannot find fission fragment yields file: " << yieldsFilename << endl; exit(-1);}

    // Set the normalization, min./max. values, and read
    double norm    = 0.0;
    double avTKE__ = 0.0;
    TKEmin = 999; TKEmax = 0; Amin = 300; Amax = 0;
    while (!yieldsFile.eof()) {
	yieldsFile >> A >> TKE >> yields;
	YATKE[A][TKE] = yields;
	norm += yields;
	if (yields != 0.0) { // Update the min./max. values if the yield is nonzero
	    if (TKE<TKEmin) TKEmin = TKE;
	    if (TKE>TKEmax) TKEmax = TKE;
	    if (A<Amin) Amin = A;
	    if (A>Amax) Amax = A;
	}
    }

    yieldsFile.close();

    Amin = Amin - 5; // To be safe
    Amax = Amax + 5;
    norm = 0.0; // TODO: I THINK WE NEED TO USE THE NORMALIZATION FROM THE SYMMETRIZED YIELDS (P.J.)
    // Symmeterize the yields
    for (int i=Asym; i<=Amax; i++) {
	for (int j=0; j<NUMTKE; j++) {
	    YATKE[i][j] = YATKE[Ac - i][j];
	    norm += YATKE[i][j];
	}
    }

    // Renormalize yields to 1
    for (int i=0; i<NUMA; i++) {
	for (int j=0; j<NUMTKE; j++) {
	    YATKE[i][j] /= norm;
	}
	for(int j=0; j<NUMTKE - 1; j++) avTKE__ += 0.5*(YATKE[i][j] + YATKE[i][j + 1])*(j + 0.5); // Integrate
    }
    cout << "Average TKE is " << avTKE__/2. << endl;
    return;
}

/*******************************************************************************
 * readYieldsAZTKE
 *------------------------------------------------------------------------------
 * Reads fission fragment yields in the form of [A, Z, TKE, Yield].
 ******************************************************************************/
void FissionFragments::readYieldsAZTKE (string yieldsFilename) {

    ifstream yieldsFile;
    int A, Z, TKE;
    double yields;
    char c[300];
    string line;

    yieldsFile.open(&yieldsFilename[0]);
    if (!yieldsFile) {cerr << "ERROR: Cannot find fission fragment yields [A,Z,TKE,Y] file: " << yieldsFilename << endl; exit(-1);}

    double norm=0.0;
    TKEmin=999; TKEmax=0; Amin=300; Amax=0; Zmin=100; Zmax=0;

    yieldsFile.getline(c,299);
    while (!yieldsFile.eof()) { // Loop through the yields file and input the necessary info
		if (c[0] == '#') {yieldsFile.getline(c,299); continue;} // Ignore comment lines

		line = (string)c;
		std::istringstream iss(line);
		int Nparam = std::distance(std::istream_iterator<string>(iss),std::istream_iterator<string>());
		if (Nparam!=4) {yieldsFile.getline(c,299); continue;} // Ignore any lines without the correct number of parameters

		// Determine the Z and A for this data
		iss.str(line);
		iss.clear();
		iss >> A;
		iss >> Z;
		iss >> TKE;
		iss >> yields;

		YATKE[A][TKE] += yields;
		YZATKE[Z][A][TKE] += yields;
		YZA2[Z][A] += yields;
		norm += yields;

		if (yields!=0.0) {
		    if (TKE<TKEmin) TKEmin=TKE;
		    if (TKE>TKEmax) TKEmax=TKE;
		    if (A<Amin) Amin=A;
		    if (A>Amax) Amax=A;
		    if (Z<Zmin) Zmin=Z;
		    if (Z>Zmax) Zmax=Z;
		}

		yieldsFile.getline(c,299);
    }

    yieldsFile.close();
    // norm = 0.0; // TODO: I THINK WE NEED TO USE THE NORMALIZATION FROM THE SYMMETRIZED YIELDS (P.J.)
    // symmetrize yields, if needed
    for (int i=Asym; i<Amax; i++) {
		for (int j=0; j<NUMTKE; j++) {
		    YATKE[i][j] = YATKE[Ac-i][j];
		    // norm += YATKE[i][j];
		}
    }
    double NewNorm = 0.0;
    //-- renormalize yields Y[A][TKE] to 1.0
    for (int i=0; i<NUMA; i++) {
		for (int j=0; j<NUMTKE; j++) {
		    YATKE[i][j] /= norm;
		    NewNorm += YATKE[i][j];
		}
    }
    // cout<<NewNorm<<endl;
    double Ymax;
    int jmax;
    // determine YZA from YZA2
    for (int i=Amin; i<=Amax; i++) {
		Ymax=-100;
		jmax=0;
		for (int j=1; j<NUMZ; j++) {
		    if (YZA2[j][i]>Ymax) {Z0[i]=j; jmax=j; Ymax=YZA2[j][i];}
		}
		Z0[i]=jmax;
    }

    for (int i=Amin; i<Amax; i++) {
		for (int j=0; j<=2*dZ; j++) {
		    YZA[j][i] = YZA2[int(Z0[i])+j-dZ][i];
		}
    }
    return;
}

/*******************************************************************************
 * spinParameter2
 *------------------------------------------------------------------------------
 * Returns the spin parameter B, which enters the spin distribution P(J). B will
 * depend on A, Z, and the excitation energy U through the moment of inertia
 * B = sqrt[I*T*hbar^2]
 * Moment of Inertia I of a rigid-body ellipsoid is given by:
 *
 * I = (2/5) * M * R^2 * [1 + 0.31*beta2 + 0.44*beta2^2 + ...]
 * M = A * m_nucleon		-- m_nucleon = (m_p + m_n)/2
 * R = r0 * A^(1/3)		-- r0 = 1.2 fm
 ******************************************************************************/
double FissionFragments::spinParameter2(double U, int A, int Z) {

    int zaid = 1000*Z+A;
    double b2 = getDeformation(zaid);

    // Rigid-body moment of inertia
    double momInertia = 0.4*pow((double)A,5./3.)*1.2*1.2*.5*(NEUTRONMASS+PROTONMASS)*amuMeV;
    momInertia *= (1.0 + 0.31*beta2[Z*1000+A] + 0.44*pow(beta2[Z*1000+A],2.));

    // Spin-scaling factor alphaI
    momInertia *= alphaI;

    double betaTemp = getTemperature(zaid,U);
    double Bspin2 = momInertia*betaTemp/hbarc/hbarc;

    return (Bspin2);
}

/*******************************************************************************
 * findEnergyIndex
 *------------------------------------------------------------------------------
 * Finds the index value k0 that corresponds to a value x0 in the array xarray
 ******************************************************************************/
template <int N> int findEnergyIndex(double x0, double (&xarray)[N]) {

    int k0 = -1;
    int numberElements = N;
    for (int k=1; k<numberElements; k++) {
	if (xarray[k]>x0) {
	    k0 = k - 1;
	    break;
	}
    }
    return k0;
}

/*******************************************************************************
 * constructYields
 *------------------------------------------------------------------------------
 * Constructs the yields Y(A), Y(Z), Y(A,Z), <TKE>(A) for a given fission event
 * after sampling whether or not a pre-fission neutron is emitted.
 ******************************************************************************/
void FissionFragments::constructYields(Yields *ffy) {

    // Initialize to the original compound system values -- these are the ones provided at the command line
    incidentEnergy = incidentEnergy0;
    Ac = Ac0;
    Zc = Zc0;
    SnCompound = SnCompound0;

    double excitEnergy = ncl[0].max_energy; // Maximum energy that is available for neutron emission
    getPrefissionNeutronEnergy(&excitEnergy); 

    // For neutron-induced fission, we shift to the new compound system
    if (!sf_flag) {
	Ac = Ac0 - n_pfn;
	massExcessCompound = mass_excess(Zc, Ac);
	SnCompound = -mass_excess(Zc,Ac) + mass_excess(Zc,Ac - 1) + mass_excess(0,1);
	incidentEnergy = excitEnergy - SnCompound; // Calculate the "effective" incident neutron energy
    }

    double WahlPE = incidentEnergy0 + SnCompound0;

    // Set the spin strength (i.e. alphaI) for this new compound nucleus and "effective" incident neutron energy
    setSpinStrength(ZAIDc - n_pfn,incidentEnergy);
//    setSpinStrength(incidentEnergy); // IS parameterization

    // Set the fissioning system for this new compound nucleus and "effective" incident neutron energy
    ffy->setFissioningSystem(ZAIDc - n_pfn,incidentEnergy);

    // Determine the new Amin and Amax
    Amin = 300; Amax = -1;
    maxYieldA = 0.;
    for (int a=0; a<NUMA; a++) {
	YA[a] = ffy->yieldA(a);
	if (YA[a]>0. && a<Amin) Amin = a;
	if (YA[a]>0. && a>Amax) Amax = a;
    }

    // Renormalize the mass distribution
    double s = 0.;
    for (int a=Amin; a<=Amax; a++){
	s += YA[a];
    }
    for (int a=Amin; a<=Amax; a++){
	YA[a] /= s;
	if (YA[a]>maxYieldA) maxYieldA = YA[a]; // Find the new maxY(A)
    }

    // Rescale the <TKE> for this new compound nucleus and "effective" incident neutron energy
    ffy->setRescaleTKE(YA,Amin,Amax);
    // Initialize and compute the Wahl parameters for a Z-distribution
    ZDistribution WahlParams = ZDistribution(Z0,sZ0,FZZ0,FNZ0);
    if (WahlParams.ReturnMethod() == "LegacyCGMF") {
    	WahlParams.ComputeZDistribution(Zc,Ac,incidentEnergy+SnCompound,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    } else {
    	WahlParams.ComputeZDistribution(Zc0,Ac0,incidentEnergy0+SnCompound0,Amin,Amax,dZ,Z0,sZ0,FZZ0,FNZ0,&Zmin,&Zmax);
    }
    // computeWahlParametersForChargeDistribution(); // Compute the new Wahl parameters for the charge distribution
    buildYZA(); // Build the Y(Z|A) and Y(Z)
    buildYZ();

    // Reset the maxY(Z)
    for (int i=0; i<NUMA; i++) {
	maxYieldZ[i] = 0.0;
	for (int j=0; j<NUMZ; j++) {
	    if (YZA2[j][i]>maxYieldZ[i]) {maxYieldZ[i] = YZA2[j][i];}
	}
    }

    return;
}





// void FissionFragments::readMultichanceFissionData(void) {

//   ifstream fp;
//   string str = MULTICHANCEFISSION;
//   string line;

//   int ZAID, ne, nn;

//   // Try current directory
//   fp.open(&str[0]);
//   // Then system data area
//   if (!fp) {
//     str = datadir + str;
//     fp.open(&str[0]);
//   }
//   if (!fp) cgmTerminateCode ("Multi-chance fission data file not found\n");

//   int numZAID = 0;
//   while (getline(fp,line)) {
//     if (line[0]=="#" or line=="") continue;
//     isstringstream(line) >> ZAID >> ne >> nn;




//     numZAID++;
//     if (numZAID>MAX_ZAIDc) cgmTerminateCode("Increase size of multi-chance fission array\n");
//   }

//   return;

// }



/*******************************************************************************
 * readMultichanceFissionData
 *------------------------------------------------------------------------------
 * Reads the multi-chance fission probability data from the input file and
 * stores the information. Creates the barrier[] and emissprob[] arrays.
 ******************************************************************************/
void FissionFragments::readMultichanceFissionData(void) {

    // TODO: RECODE THIS SO IT DOESN'T HANG IF THE MULTICHANCEFISSION.DAT FILE HAS THE INCORRECT FORMAT -- MIGHT NEED TO RESTRUCTURE THE FILE TOO
    // -----------------------------------------------------------------------------------------------------------------------------
    string fileFissionData = datadir;
    fileFissionData += "multichancefission.dat";
    ifstream fissdata(fileFissionData.c_str(),ios::in);

    if (!fissdata.is_open()) {
	cerr << "[readMultichanceFissionData] file " << fileFissionData << " was not found" << endl; exit(-1);
    }

    char c[400];
    string line, subline;
    fissdata.getline(c,399);
    bool found = false; // Tags when we've determined the multi-chance fission probabilities for this incident particle energy

    while (!fissdata.eof()) {

	int id, nenergies, nn;

	if (c[0] == '#') {fissdata.getline(c,399); continue;} // Ignore comment lines
	line = (string)c;
	std::istringstream iss(line);

	int Nparam = std::distance(std::istream_iterator<string>(iss),std::istream_iterator<string>());
	if (Nparam<3) {fissdata.getline(c,399); continue;} // Ignore empty lines (or those with less than 3 items)

	iss.str(line); // Reload the line
	iss.clear();
	// Load the ZAID, number of lines, and max number of neutrons to emit
	iss >> id; iss >> nenergies; iss >> nn;

	if (ZAIDc - 1 == id) { // Neutron-induced fission
	    barrier = new double[nn + 1]; // Array for barrier heights
	    for (int i=0;i<=nn;i++) {
		fissdata >> barrier[i]; // Read in barrier heights
	    }
	    emissprob = new double[nn + 1]; // Array for final emission probabilities (at En = incidentEnergy)
	    double Pf[nenergies][nn+1]; // Fission probabilities for each energy and multi-chance
	    double energies[nenergies]; // Value of the incident particle energy

	    // Read all the probabilities and energies
	    for (int i=0;i<nenergies;i++) {
		fissdata >> energies[i];
		for (int j=0;j<=nn;j++) fissdata >> Pf[i][j];
	    }

	    if (incidentEnergy<energies[0] || incidentEnergy>energies[nenergies-1]) {
		cerr << "[readMultichanceFissionData] could not find multi-chance fission data for this incident energy" << endl;
		exit(-1);
	    }

	    // Find the multi-chance fission probabilities at the specified incident particle energy
	    for (int i=0;i<nenergies;i++) {
		if (incidentEnergy<=energies[i]) {
		    double slope;
		    for (int j=0;j<=nn;j++) {
			slope = (Pf[i][j] - Pf[i-1][j])/(energies[i] - energies[i-1]);
			if (slope!=0.0) {
			    emissprob[j] = (incidentEnergy - energies[i-1])*slope + Pf[i-1][j];
			} else {
			    emissprob[j] = Pf[i][j];
			}
		    }
		    break;
		}

	    }

      // PT, April 2021: segmentation fault if nemit>=nn+1... simple FIX for now: replace <=nemit with <=min(nemit,nn)

	    // Determine the maximum number of pre-fission neutrons
	    int k = nn;
//      for (int j=0;j<=nemit;j++) { // TODO: If "statTotalNeutrons" returns a value lower than that in the multi-chance file, we'll only use the "statTotalNeutrons" value.
      for (int j=0;j<=min(nemit,nn);j++) { // TODO: If "statTotalNeutrons" returns a value lower than that in the multi-chance file, we'll only use the "statTotalNeutrons" value.
		if (emissprob[j]<=0.0) {
		    k = j - 1;
		    break;
		}
	    }
	    nemit = k;

	    // Determine the maximum probability (for sampling later)
	    emissprob_max = 0.0;
	    double s = 0.0; // Normalize
	    for (int j=0;j<=nemit;j++) {
		s += emissprob[j];
	    }
	    for (int j=0;j<=nemit;j++) {
		emissprob[j] /= s;
		if (emissprob[j]>emissprob_max) emissprob_max = emissprob[j];
	    }
	    found = true;
	    break;

	} else { // Move to the next nucleus
	    for (int k=0;k<=nenergies + 1;k++) {
		fissdata.getline(c,399);
	    }
	}
    }

    fissdata.close();

    // If we didn't find the multi-chance fission probabilities and we need it, then exit
    if (!found && nemit > 0) {cerr<<"[readMultichanceFissionData] No multichance fission probabilities for this nucleus!"<<endl; exit(-1);}

    return;
}

/*******************************************************************************
 * getNumberPrefissionNeutrons
 *------------------------------------------------------------------------------
 * Given an excitation energy, sample whether we emit a pre-fission neutron or not 
 ******************************************************************************/
int FissionFragments::getNumberPrefissionNeutrons(void) {


   int n;
   do { // Sample a number of pre-fission neutrons to emit
	n = (int) floor(rng_cgm()*(nemit + 1));
   } while (emissprob[n] < emissprob_max*rng_cgm());

   return n;

}

/*******************************************************************************
 * getPrefissionNeutronEnergy
 *------------------------------------------------------------------------------
 * Given an excitation energy, sample whether we emit a pre-fission neutron or
 * not and determine its energy. For each pre-fission neutron emitted, we
 * calculate their energy (energy_pfn), spectrum (spn_pfn), and multiplicity-
 * dependent spectrum (spn_pfn_mult).
 ******************************************************************************/
void FissionFragments::getPrefissionNeutronEnergy(double *excitEnergy){

    // If the max. number of neutrons to emit is 0, then move on
    if (nemit == 0) {
	return;
    }

    int n = n_pfn;

    //int n;
    //do { // Sample a number of pre-fission neutrons to emit
//	n = (int) floor(rng_cgm()*(nemit + 1));
    //} while (emissprob[n] < emissprob_max*rng_cgm());

    // If we sampled 0 pre-fission neutrons, move on
    if (n == 0) return;

    // Determine the maximum pre-fission neutron energy (in c.o.m. frame)
    double eps_max = ncl[n].excitation[0] - barrier[n];
    // Can't sample negative energies -- i.e. we're below the barrier of the 2nd-chance daughter
    if (eps_max < 0.) preeq_flag = false;
    if (eps_max < 0.) return;

    int kfmax = (int) ceil(eps_max/ENERGY_BIN);

    // Sample the first pre-fission neutron from the spectrum sp1.dat (either evaporation or pre-fission)
    double max_prob = 0.;
    for (int i=0; i<kfmax; i++){
	if (sp1[i] > max_prob) max_prob = sp1[i]; // Update the maximum probability for the spectrum sp1
    }

    // Sample the first neutron-out spectrum 'sp1' (contains evaporation or pre-equilibrium contributions)
    int kf;
    do { // Sample from the spectrum
	kf = (int) floor(rng_cgm()*kfmax);
    } while (sp1[kf] < max_prob*rng_cgm());

    cmEnergy_pfn[0] = energy_grid1[kf] + 0.5*ENERGY_BIN*(1-2*rng_cgm()); // Continuum-to-continuum energy transfer

    // if the neutron was pre-equilibrium, construct the CDF for the angular sampling
    float Eexc = ncl[1].excitation[kf];
    float Einc = incidentEnergy;
    if (preeq_flag){
       float p0 = 0., p1 =0.;
       float aparm = 0., bparm = 0., cparm = 0., dparm = 0.;
       for (int i=0;i<6;i++){
	  float Epow = pow(Eexc,i);
	  aparm = aparm + preeqScatteringParams[0][0][i]*Epow;
	  bparm = bparm + preeqScatteringParams[0][1][i]*Epow;
	  cparm = cparm + preeqScatteringParams[0][2][i]*Epow;
	  dparm = dparm + preeqScatteringParams[0][3][i]*Epow;
       } 
       p0 = aparm + bparm*Einc + cparm*Einc*Einc + dparm*Einc*Einc*Einc;
       aparm = 0.;
       bparm = 0.;
       cparm = 0.;
       dparm = 0.;
       for (int i=0;i<6;i++){
	  float Epow = pow(Eexc,i);
	  aparm = aparm + preeqScatteringParams[1][0][i]*Epow;
	  bparm = bparm + preeqScatteringParams[1][1][i]*Epow;
	  cparm = cparm + preeqScatteringParams[1][2][i]*Epow;
	  dparm = dparm + preeqScatteringParams[1][3][i]*Epow;
       } 
       p1 = aparm + bparm*Einc + cparm*Einc*Einc + dparm*Einc*Einc*Einc;

       // construct the angular distribution - p0*sinh(cos(theta))+p1
       // then construct the CDF
       double cosTheta;
       Ppreeq[0][0] = -1.0;
       Ppreeq[0][1] = 0.0;
       for (int i=1;i<THETA_STEPS;i++){
	  cosTheta = 2.0*i/(THETA_STEPS-1)-1.0;
          if (p0==0. && p1==0.) {
             Ppreeq[i][0] = cosTheta;
             Ppreeq[i][1] = 1.;
          } else {
	     Ppreeq[i][0] = cosTheta;
	     Ppreeq[i][1] = Ppreeq[i-1][1] + p0*sinh(cosTheta)+p1;
          }
       }
       // normalize the distribution so it is a CDF
       for (int i=0;i<THETA_STEPS;i++){
	  Ppreeq[i][1] = Ppreeq[i][1]/Ppreeq[THETA_STEPS-1][1];
       }
    }

    // Sample higher-order multi-chance fission probabilities (evaporation spectra only)
    int kstart = kf;
    for (int i=2;i<=n;++i) {

	evapInterface(pdt+neutron,tc,i-1,kstart,sp);
	eps_max = ncl[n].excitation[kstart] - barrier[n];
	kfmax = kstart + (int)floor(eps_max/ENERGY_BIN) + 1;

	if (eps_max<0.) { // << 1.0.6 >> Can this happen? ...
	    cout << i << " nemit_now=" << n << " *eps_max=" << eps_max << " kstart=" << kstart << endl;
	}

	if( eps_max<ENERGY_BIN){ // IS: added a trap to catch zero energy prefission neutrons
	  cmEnergy_pfn[i-1]=1e-10; // IS: changed to avoid zero-energy neutrons
	  continue;
	}

	max_prob = 0.;
	for (int k=kstart;k<kfmax;k++) {
	    if (sp[k]>max_prob) max_prob = sp[k];
	}

	do {
	    kf = (int) (kstart + (kfmax - kstart)*rng_cgm());
	} while (sp[kf]<max_prob*rng_cgm());

	if (kf<0) cout << i << " nemit_now=" << n << " kstart=" << kstart << " kf=" << kf << endl;

	cmEnergy_pfn[i-1] = tc[kf].ecms + 0.5*ENERGY_BIN*(1-2*rng_cgm());
	kstart = kf;

    }
    *excitEnergy = ncl[n].excitation[kf];

    return;
}

/*******************************************************************************
 * getSpinScalingCoefficients
 *------------------------------------------------------------------------------
 * Given a compound nucleus ZAID and the spontaneous fission flag, locate the
 * spin-scaling coefficients alpha_0 and alpha_slope. If none are found, then
 * leave the input values unchanged.
 ******************************************************************************/
void FissionFragments::getSpinScalingCoefficients(int ZAID_0, bool sf_flag, double *alpha_0, double *alpha_slope) {

	int zaid = ZAID_0;
	// For (sf) reactions, the ZAID is -ZAID
	if (sf_flag) {
		zaid = - zaid;
	}
	// Loop through the alpha(En) coefficients
	for (int i=0;i<ALPHA_MAX_ZAIDc;i++) {
		// Check if the ZAID in the table matches this reactions ZAID
		int table_zaid = int(alpha_coefficients[i][0]);
		if (zaid == table_zaid) {
			// Set the coefficients
			*alpha_0 = alpha_coefficients[i][1];
			*alpha_slope = alpha_coefficients[i][2];
			break;
		}
	}

	return;
}

/*******************************************************************************
 * setSpinStrength
 *------------------------------------------------------------------------------
 * Sets the spin-scaling factor alphaI for a specified "effective" incident
 * neutron energy einc and the compound nucleus ZAID_0. TODO: PUT THIS INFORMATION IN A DATA FILE (P.J.)
 ******************************************************************************/
void FissionFragments::setSpinStrength(int ZAID_0, double einc) {

    // Default values -- from systematics of known systems
    double alpha_slope = 0.071;
    double alpha_0 = 1.5;
    // Spontaneous fission reactions have a different set of defaults
    if (sf_flag) {
    	alpha_slope = 0.0;
    	alpha_0 = 1.3;
    }

    if (alphaI0>0.) {
		alphaI = alphaI0; // Always use the user-defined alphaI if it was provided
		return;
    }

    // Assign the spin-scaling coefficients
    getSpinScalingCoefficients(ZAID_0,sf_flag,&alpha_0,&alpha_slope);

    if (einc<0.) {
		alphaI = alpha_0;
    } else {
		alphaI = alpha_0 + alpha_slope*einc;
    }

    return;
}

// IS: updated the spin dependence on the excitation energy
void FissionFragments::setSpinStrength(double einc){

  switch(ZAIDc0){
  case 94240:


    //    alphaI=2.57-1.07*exp(-0.15*einc);

    if(incidentEnergy0==0.0)
      alphaI=1.5;
    else
      alphaI=2.16-.633*exp(-.233*einc); // old parameterization

    break;

  case 94242:

    if(incidentEnergy0==0.) 
      alphaI=1.5;

    break;

  case 92236:

    if(einc>=1.)
      alphaI=2.3-1.474*exp(-0.47*einc);
    else
      alphaI=1.624-0.25*einc;

    break;

  case 92239:

    /* Original

    if(einc<1.5)   
      alphaI=1.44;
    else
      alphaI=2.49-7.11*exp(-1.28*einc); // this would be negative at around 0.8, keep alphaI constant below 1.5 MeV
    break;
    */

    //    alphaI=3.42-7.52*exp(-1.28*einc);
    //    alphaI=3.91-57.72*exp(-2.015*einc);
    break;
    if(einc>=2.)
      alphaI=3.8-40.58*exp(-1.8286*einc);
    else
      alphaI=3.8-40.58*exp(-1.8286*2.);
    break;

  default:
    break;
  }

}

void readPreeqData(void){

  ifstream fp;
  string str = PREEQDATAFILE;

  // Try current directory
  fp.open(str.c_str());

  // Then system data area
  if (!fp) {
    str = datadir + str;
    fp.open(str.c_str());
  }

  if (!fp) cgmTerminateCode("preequilibrium data file not found");

  string line;

  while(true){
    getline(fp,line);
    string c=line.substr(0,1);
    if(c!="#"){
      n_reactions = atoi(c.c_str());
      break;
    }
  }

  for(int i=0; i<n_reactions;i++)
    fp >> preeq_id[i] >> preeq_params[i][0] >> preeq_params[i][1] >> preeq_params[i][2] >>preeq_params[i][3] ;

  fp.close();

  return;

}





