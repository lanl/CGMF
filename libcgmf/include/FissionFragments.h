/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file FissionFragments.h

  \brief Scission fragment yields, either read in or from systematics 

*/

#ifndef __FISSIONFRAGMENTS_H__
#define __FISSIONFRAGMENTS_H__

#include <iostream>
#include <string>
#include <vector>

#include "config-ff.h"
#include "Yields.h"
#include "structur.h"

using namespace std;

// Data file names
#define MASSESFILE        "masses.dat"
#define DEFORMATIONSFILE  "deformations.dat"
#define TEMPERATURESFILE  "temperatures.dat"
#define LEVELDENSITYFILE  "ldp.dat"
#define MULTICHANCEFILE   "multichancefission.dat"
#define ANISOTROPYFILE    "anisotropy.dat"
#define PREEQSCATTERFILE  "preeq-scattering-params.dat"
#define PREEQDATAFILE     "preeq-params.dat"
#define SPINSCALINGFILE   "spinscalingmodel.dat"
#define TKEYIELDFILE      "tkemodel.dat"
#define MASSYIELDFILE     "yamodel.dat"
#define RTAFILE           "rta.dat"

// Accessors 
void   readMasses         (void);
double getMassExcess      (int);
void   readDeformations   (void);
double getDeformation     (int);
void   readTemperatures   (void);
double getTemperature     (int, double);
void   readLDP            (void);
double getLDP             (int, double);
void   readAnisotropy     (void);
void   readPreEqAngularDistributionParameters (void);
void   readPreeqData      (void);
void   readSpinScaling    (void);
void   readRTA            (void);

#define MAX_ZAID 140000
#define NUM_GRID 32
static double masses                [MAX_ZAID];
static double beta2                 [MAX_ZAID];
static double temperatures          [MAX_ZAID][NUM_GRID];
static double temperatureEnergyGrid [NUM_GRID];
static double ldp2                  [MAX_ZAID][NUM_GRID];

// Anisotropy
#define MAX_ZAIDc 7
#define NUM_EGRID 201

static double anisotropy[MAX_ZAIDc+1][NUM_EGRID+1] = {{1.0}};
static double Pcostheta[THETA_STEPS][2];

// Spin scaling (alpha) dependence on En
#define ALPHA_MAX_ZAIDc 15
static double alpha_coefficients[ALPHA_MAX_ZAIDc][3];

// Multi-chance fission probabilities
//static double multiChanceFissionProbabilities [MAX_ZAIDc][5][41]; // ZAID, 1st, 2nd, 3rd, 4h chance probabilities as a function of energy

// Preequilibrium fit parameters for inelastic scattering distributions
static double preeqScatteringParams[2][4][6];
static double Ppreeq[THETA_STEPS][2];

//Preequilibrium data for up to MAX_TARGETS reactions
#define MAX_TARGETS 10  //maximum number of reactions handled
static int n_reactions; // the number of reactions in preeq-params.dat file
static int    preeq_id    [MAX_TARGETS];
static double preeq_params[MAX_TARGETS][4];

template<int N> int findEnergyIndex (double x0, double (&xarray)[N]);

struct fissionFragmentType {
  int Z, A, N;
  double mass; // ground-state mass (MeV)
  double U;
  double spin;
  int parity;
  double KE;
  double KEpost;
  double TKE, energyRelease; // also depend on complementary fragment
	double preMomentum[3]; // (classical) momentum components (x,y,z) in LAB frame
	double postMomentum[3]; // after neutron emission
};

struct RTAtype {
	int ZAID;
	double ratios [NUMA];
	RTAtype () {
		std::fill_n(ratios, NUMA, 1.0);
	}
};

static int NUMZAID_RTA;
static RTAtype *RTAdata;
void cleanupRTA (void);


class FissionFragments
{
  
  /* PUBLIC */
  public :
	
	// FOR TESTING TEMPERATURES PURPOSE ONLY ---------------------------------------------------- /////////////////////
	void computeTemperatures (int, int);
	void readTemperatures (void);
//	void readLevelDensityParameters (void);
	
	double temp [2000][32];
	double ldp [2000][32];
	int ZAIDlist [2000];
	
  FissionFragments (int ZAID, double incidentEnergy, double alphaSpin);
  FissionFragments (int ZAID, double incidentEnergy, double *alphaSpin, string yieldsFile);
  FissionFragments (int ZAIDf, double excitationEnergy, const string);
  FissionFragments (void);
	
  ~FissionFragments (void);
  
  void writeYieldsInFile (string outputFilename);
  void generateInitialFissionFragmentHistories (string outputFilename, int numberEvents);
  void generateInitialFissionFragmentHistories (fissionFragmentType*, fissionFragmentType*, const int numberEvents);
  void generateSingleFissionFragments (fissionFragmentType*, fissionFragmentType*, int ZAIDf, double Eexc, int numberEvents);
  void checkDistributions (string inputFilename, string outputFilename);
  
  void studyEnergySorting (void);

  void sortExcitationEnergy (double Utot, int Al, int Ah, double *Ul, double *Uh, double RT); // for Carjan

  int Amin, Amax;       // min and max of fragment mass
  int Zmin, Zmax;       // min and max of fragment charge
  int TKEmin, TKEmax;   // min and max of total kinetic energy (MeV)
  int Zpmin, Zpmax;     // min and max of most probable charges
  
  int ZAIDc;            // ZAID of compound (fissioning) nucleus
  int ZAIDt;            // ZAID of target nucleus
  int Ac, Zc;           // mass and charge of compound (fissioning) nucleus
  int At, Zt;           // mass and charge of target nucleus
  int Ap, Zp;           // mass and charge of projectile
  int Asym, Zsym;       // mass and charge of symmetric fragment
	
	int ZAIDc0, Ac0, Zc0; // Original (before any pre-fission neutron emission) compound fissioning nucleus

  double incidentEnergy;  // incident neutron energy [MeV] or equivalent incident neutron energy if a neutron is emitted
	double incidentEnergy0; // incident neutron energy [MeV] BEFORE any pre-fission neutron emission
	
  double YA[NUMA];
  double YATKE[NUMA][NUMTKE];
  double YTKE[NUMTKE];
  double YZ[NUMZ];
  double YZA[NUMdZ][NUMA];
  double YZA2[NUMZ][NUMA];

  double YTKEf[NUMTKE];

  // Average Q-values as a function of TKE and A
  double QfTKE[NUMTKE];
  double QfA[NUMA];
  
  double maxYieldZ[NUMA];
  double maxYieldATKE;
	
	double maxYieldA;	// << 1.0.6 >>
	
  double YZATKE[NUMZ][NUMA][NUMTKE];

  // AL -- Anisotropy CDF
  void computeFragmentAngularDistribution (void);
  double anisotropyCoefficient;
	
	
  /* PRIVATE */
  private :
  
  void init (void);
  void setOptions (void);
	
  void readMainInputFile (string inputFilename);
//  void createDefaultInputFile (string inputFilename, int ZAID, double incidentEnergy, double alphaSpin);
  
  // To include with XML I/O
  //  void readXMLInputFile (string inputFilename);
  //	void createDefaultXMLInputFile (string inputFilename, int ZAID, double incidentEnergy);
  
  void buildYields (void);
  void readYieldsATKE (string);
	void readYieldsAZTKE (string);
  
  void sampleFissionFragments (fissionFragmentType *lf, fissionFragmentType *hf);

	// << 1.0.6 >>
	void sampleFissionFragments (fissionFragmentType *lf, fissionFragmentType *hf, int);
	void sampleFissionFragments (fissionFragmentType *lightFragment, fissionFragmentType *heavyFragment, Yields *ffy) ;
	
	// << 1.0.7 >> added fragment momentum calculation
	void computeFragmentInitialConditions (double TKE, int Zl, int Al, int Zh, int Ah, double *Ul, double *Uh,
																				 double *spinl, double *spinh, int *parl, int *parh, double *energyRelease,
																				 double *momentum_l, double *momentum_h);
	
  void computeFragmentExcitationEnergy (double TXE, int Zl, int Al, int Zh,
                                        int Ah, double *Ul, double *Uh, string sortingOption);
  
  void getInitialSpin (int Z, int A, double *U, double *spin);
  
  void computeWahlParametersForChargeDistribution (void);
  void buildWahlYA (void);
  void buildYATKEfromSystematics (void);
  void buildSystematicsTKEvsMass (void);
	
	// << 1.0.6 >>
  void buildYATKEfromSystematics (Yields *);
  void buildYATKE (Yields *);
	
  void buildYA (void);
  void buildYZA (void);
  void buildYZ (void);
  void buildYTKE (void);
  
  void resetYields (void);
  
  // energy sorting methods & variables
  double Pfactor (int Z, int N);
  void computePfactors (void);
#ifdef DEVUTIL
  void computeLevelDensityTables (void);
  void computeLevelDensityParameterTables (void);
#endif

  void readLevelDensityParameterTables (string filename);
	void readLevelDensityTables (string filename);
    
  double computeEnergyFromMaxEntropy (int Zl, int Al, double TIXE);
  
  void readRTAh (int Zt, int At);
  void readRTAh (int Zt, int At, string filename);
  void readMollerDeformationEnergies (string filename);  
  
  double Pfactors[NUMA][NUMZ];
  double RTAh[NUMA];
  
  // level density tables f(A,Z,E)
  double ***levelDensities;
  double ***levelDensityParameters;
  double ldEnergyGrid [NUME];
	double ldEnergyGrid2 [NUME2];
  
  //-- IONEL
  //	double beta2 [NUMZ][NUMA]; // put in PUBLIC 
  void readDeformations (void);
  
  double ***temperatures;
  double spinParameter2 (double, int, int);
  void setSpinStrength (double); // Ionel, new version of high-Einc yields
  void setSpinStrength (int,double); // PJ, parameterization similar to <TKE>
  void getSpinScalingCoefficients(int,bool,double*,double*);
  void computeTemperatureTables (void);
  
//  Nucleus *fragments; // light and heavy fragments following CGM structure
  
  void readNeutronBindingEnergies (void);
//  void readMasses (void);
  void computeNeutronSeparationEnergies (void);
  
  // Wahl's parameters
  double Z0[NUMA];
  double sZ0[NUMA];
  double FZZ0[NUMA];
  double FNZ0[NUMA];
  
  double sigmaZ;
  int dZ, dZ1;
  
  double meanTKE[NUMA];
  double sigTKE[NUMA];
  
  string energySortingOption;
  double RT;
	double ratioU;
  double alphaI0,alphaI;

  double EdefA [NUMA]; // deformation energies (MeV)
  
  int lightFragmentSpinParameter;
  int heavyFragmentSpinParameter;
  
  string yieldsOption;
  string yieldsFilename;
  
  double B1n[NUMZ][NUMA];
  double B2n[NUMZ][NUMA];
  
  double S1n[NUMZ][NUMA];
  double S2n[NUMZ][NUMA];
  
  double ExpMassMeV[NUMZ][NUMA]; // excess masses in amu
  double ExpMassAMU[NUMZ][NUMA]; // excess masses in MeV
  
  double massExcessCompound; // excess mass in MeV for the compound fissioning nucleus
  double SnCompound; // neutron separation energy for the compound fissioning nucleus

	// << 1.0.6 >>
	
	double SnCompound0; // neutron separation energy for the original (before any pre-fission neutron emission) compound fissioning nucleus
	
	void constructYields (Yields *ffy);        // constructs yields for multi-chance fission
	void cgmBinaryNeutron (int nemit,int c0);  // this is a wrapper similar to specMain
	
	// pre-fission neutron emission << 1.0.6 >>
	
	double *sp1;   // total spectrum for the first neutron emitted
	double *sp;    // evaporation spectrum from the neutrons
	Pdata pdt[MAX_CHANNEL];
	Transmission *tc;
	int nemit;
	double * barrier, * emissprob, emissprob_max;

	void readMultichanceFissionData (void); // void read_multifission_data(void);
	void getPrefissionNeutronEnergy (double *); // int get_prefission_neutrons(double *);
	int getNumberPrefissionNeutrons (void);

	double *energy_grid1;
	double *sp_pfn;
  double cmEnergy_pfn[MAX_NUMBER_NEUTRONS_PRE] ;// the energy of the generated prefission neutrons
  double labEnergy_pfn[MAX_NUMBER_NEUTRONS_PRE]; // energy of the boosted pre-fission neutrons
  int neutronType_pre[MAX_NUMBER_NEUTRONS_PRE]; // this is a placeholder, to be used later to identify
                                                // either pre-equilibrium or statistical neutrons
  double cmVelocities_pfn [MAX_NUMBER_NEUTRONS_PRE][3];
  double preeq_fraction(void);  // fit to CoH 3.5.3 pree-quilibrium fraction

  bool sf_flag; // Flag for spontaneous fission
  bool anisotropy_flag; // flag for the existance of a fit for anisotropy
  bool preeq_flag = false; // flag for sampled pre-equilibrium neutron
  int n_pfn = -1; // number of pre-fission neutrons to emit

	double prCheck[10];
	int **sp_pfn_mult;
	
	int index_save;
	
	bool recomputeYields;
	
  public:
  
  double labVelocities_pfn [MAX_NUMBER_NEUTRONS_PRE][3];

  
}; // end class Fission Fragments

/* Charge distribution class - follows the Wahl coded systematics. All massess are primary fragment masses */
class ZDistribution {

  public:

    /* Constructor and destructor */
    ZDistribution(double*,double*,double*,double*);
    ~ZDistribution();

    /* Equations 17 and 19 of Wahl's systematics */
    double WahlEq17(int,int,double,vector<double>);
    double WahlEq19(int,int,double,vector<double>);

    /* Computes the primary fragment mass boundaries */
    void ComputeMassBounds(int,int,double);

    /* Computes the Zp parameters */
    void ComputeZpParams(int,int,double);

    /* Computes the charge distribution parameters for each mass region */
    void PeakRegion(int,int,int);
    void SymmetryRegion(int,int,int);
    void WingRegion(int,int,int);
    void FarWingRegion(int,int,int);

    /* Computes the charge distribution parameters delZ(A), sigZ(A), FN(A), FZ(A) */
    void ComputeZDistribution(int,int,double,int,int,int,double*,double*,double*,double*,int*,int*);

    /* Returns the CGMF method */
    string ReturnMethod();

  private:

    /* Calculation Method */
    string Method = "NewCGMF"; // "NewCGMF" or "LegacyCGMF"
    // string Method = "LegacyCGMF"; // "NewCGMF" or "LegacyCGMF"

    /* Temporary internal Zp parameterization */
    double delz[NUMA];
    double sigz[NUMA];
    double Fz[NUMA];
    double Fn[NUMA];

    /* Compound nucleus excitation energy separating low- from mid- from high-energy fission */
    double PE_L, PE_H;

    /* Wahl parameters for the Z distribution - see Tab. 2 and 3 */
    double sigz140;
    double delz140;
    double Fz140;
    double Fn140;
    double sigzSL;
    double delzSL;
    double FzSL;
    double FnSL;
    double SL50;
    double sigz50;
    double delzmax;
    double sigzSLW;
    double delzSLW;
    double FzSLW;
    double FnSLW;

    /* Same as above but evaluated at the low- and high-energy boundaries - used in the linear interpolation at mid-energy */
    double sigz140_L, sigz140_H;
    double delz140_L, delz140_H;
    double Fz140_L  , Fz140_H;
    double Fn140_L  , Fn140_H;
    double sigzSL_L , sigzSL_H;
    double delzSL_L , delzSL_H;
    double FzSL_L   , FzSL_H;
    double FnSL_L , FnSL_H;
    double SL50_L   , SL50_H;
    double sigz50_L , sigz50_H;
    double delzmax_L, delzmax_H;
    double sigzSLW_L, sigzSLW_H;
    double delzSLW_L, delzSLW_H;
    double FzSLW_L  , FzSLW_H;
    double FnSLW_L  , FnSLW_H;

    /* Coefficients for the above parameters - see Tab. 2 (low-energy L) and 3 (high-energy H) */
    vector<double> sigz140_Lc, sigz140_Hc;
    vector<double> delz140_Lc, delz140_Hc;
    vector<double> Fz140_Lc  , Fz140_Hc;
    vector<double> Fn140_Lc  , Fn140_Hc;
    vector<double> sigzSL_Lc , sigzSL_Hc;
    vector<double> delzSL_Lc , delzSL_Hc;
    vector<double> FzSL_Lc   , FzSL_Hc;
    vector<double> FnSL_Lc   , FnSL_Hc;
    vector<double> SL50_Lc   , SL50_Hc;
    vector<double> sigz50_Lc , sigz50_Hc;
    vector<double> delzmax_Lc, delzmax_Hc;
    vector<double> sigzSLW_Lc, sigzSLW_Hc;
    vector<double> delzSLW_Lc, delzSLW_Hc;
    vector<double> FzSLW_Lc  , FzSLW_Hc;
    vector<double> FnSLW_Lc  , FnSLW_Hc;

    /* Primary (pre-neutron emission) fragment boundaries */
    double B1, B2, B3, Ba, Bb, B4, B5, B6;
    int iB1, iB2, iB3, iBa, iBb, iB4, iB5, iB6;

    /* Some values/constants at the boundaries */
    double sigzB5, delzB5, AKH;
};

#endif //__FISSIONFRAGMENTS_H__

