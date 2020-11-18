/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file specmc.cpp

  \brief Monte Carlo calculations for neutron and gamma-ray emissions

*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "mchf.h"
#include "terminate.h"
#include "config.h"
#include "ripl2levels.h"
#include "physics.h"
#include "kinematics.h"
#include "FissionFragments.h"

extern std::function< double(void) > rng_cgm;

void recordEmittedParticles (void);
void recordEmittedParticlesFission (fissionEventType *);

static void    mcHistoryLoop            (MCPoint, Pdata *, double **);
//static MCPoint mcFindNext               (int, int);
static MCPoint mcFindNext               (int, int, int, int);
static int     mcGammaCascade           (int, int, double **);
static MCPoint mcSampleStart            (double, Parity **);
static MCPoint mcFixedStart             (double, Parity **);
static void    mcRestorePopulation      (MCPoint, Nucleus *, Nucleus *);
static void    mcNormalizeSpectrum      (unsigned long, double **);
static void    mcWriteHistory           (void);

static inline void   mcStoreHistory     (bool, MCPoint, MCPoint, double **);
static inline double mcCmsToLabEnergy   (double, double);

static int       np = 0;           // number of emitted gamma-rays and neutrons
static MCHistory hs[MAX_CASCADE];  // emitted energy and particle identifier
static double    cnKE = 0.0;       // compound nucleus kinetic energy in Lab

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void    specPopCheck        (Nucleus *);
#endif

static bool firstcalcspec = true;
Transmission  tin;
Pdata         pdt[MAX_CHANNEL];
Parity        **pin = NULL;

static double time_coincidence_window=-1; // infinit time-coincidence window 

void set_time_coincidence_window(double timew){
  time_coincidence_window=timew;
}

double gammaTimeStamp [ MAX_NUMBER_GAMMAS ];

/*!
 *
 *  \brief allocate tin, pin
 *
 *  */
void specMCalloc(void)
{
    if (firstcalcspec) {
      try{
        tin.tran = new double [3*MAX_J];
        pin      = new Parity *[MAX_ENERGY_BIN];
        for(int j=0 ; j<MAX_ENERGY_BIN ; j++) pin[j] = new Parity[MAX_J];
      }
      catch(bad_alloc){
        cgmTerminateCode("memory allocation error in specMain");
      }
      firstcalcspec=false;
    }
}

/*!
 *
 *  \brief deallocate tin, pin
 *
 *  */
void specMCcleanup(void)
{
    if (!firstcalcspec) {
      delete [] tin.tran;
      for(int j=0 ; j<MAX_ENERGY_BIN ; j++) delete [] pin[j];
      delete [] pin;
      firstcalcspec=true;
    }
}


/*!

 \brief Main Monte Carlo routine for CGMF

*/
void specMCMainFission(fissionEventType *fe)
{	

  int    Z     = fe->Z;
  int    A     = fe->A;
  double U     = fe->U;
  double initJ = fe->spin;
  int    initP = fe->parity;

  specMCalloc();

  // set initial charge, mass and excitation energy of decaying fission fragment
	ncl[0].za.setZA (Z, A);
	ncl[0].max_energy = U;
	
	// how many neutrons can be emitted
	int nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
	
	/*** initialize system parameters */
	statSetupInitSystem(nemit,pdt);
	
	/*** setup discrete levels, energy bins and level densities */
	getDiscreteLevels (ncl, nemit);
	for(int i=0 ; i<=nemit ; i++){
		statSetupEnergyBin(&ncl[i]);
		statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
	}

	/*** GDR parameters */
	GDR gdr[MAX_GDR];
	statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,0.0);
	
	/*** clear spectrum array */
	for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
		for(int i=0 ; i<SPECTRA_OUTPUT ; i++) spc[i][k] = 0.0;
	}
	
	/*** start Monte Carlo */

	/*** decay start at a given discrete level */
	if(ncl[0].ncont == 0){

		/*** find nearest discrete level for given Ex */
		int nstart = statStoreLevelExcite(&ncl[0]);
  	np = 0;
		int final_state = mcGammaCascade(0,nstart,spc);
		/*** insert final state */
		int c = 0;
		MCPoint mcurr;
		mcurr.set(c,ncl[c].lev[final_state].parity,(int)ncl[c].lev[final_state].spin,final_state,ncl[c].lev[final_state].energy);
		mcurr.setDiscrete();
		hs[np++].set(gammaray, 0.0, mcurr);
		recordEmittedParticlesFission(fe);
	}

	/*** start with the continuum bin */
	else{
		/*** set up initial population */
		statSetupInitialPopulationSpec(initJ,initP,0,0.0,&tin);
		for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
			for(int j=0 ; j<MAX_J ; j++) {
				pin[k][j] = ncl[0].pop[k][j];
			}
		}

		MCPoint minit;
		double  emax[MAX_COMPOUND];
		for(int i=0 ; i<=nemit ; i++) emax[i] = ncl[i].max_energy;
		
		/*** choose initial state */
		for(int i=0 ; i<=nemit ; i++) ncl[i].max_energy = emax[i];
		
    minit = mcFixedStart(emax[0],pin);

		/*** set CN kinetic energy if provided */
		cnKE=0; // if(targE > 0.0) cnKE = targE; // -- CHECK!
		
		/*** perform simulation */
		np = 0;
		mcHistoryLoop(minit,pdt,spc);
		recordEmittedParticlesFission (fe);

	}

}

//************************************************************************************

//
// targE :: compound nucleus kinetic energy (set to 0.0 if not used?)
// initJ :: initial spin in CN
// initP :: initial parity in CN
// spinf :: spin cut-off parameter scaling factor (set to 0.0 to not change anything)


/******************************************************************************/
/* Main Monte Carlo routine for CGM (NOT fission!)                            */
/******************************************************************************/
void specMCMain(double initJ, int initP, double targE, double spinf, unsigned long nsim, double **spc)
{
  Transmission  tin;
  Pdata         pdt[MAX_CHANNEL];
  Parity        **pin = NULL;

  specMCalloc();

//---------------------------------------
//      Parameter Setting

  /*** how many neutrons can be emitted */
  int nemit = 0;
  if(INCLUDE_NEUTRON_EMISSION){
    nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
  }

  /*** initialize system parameters */
  statSetupInitSystem(nemit,pdt);

  /*** setup discrete levels, energy bins and level densities */
//	riplReadDiscreteLevels(ncl, reassign, nemit);
	getRiplDiscreteLevels (ncl, nemit);

  for(int i=0 ; i<=nemit ; i++){
//    ncl[i].ndisc = riplReadDiscreteLevels(&ncl[i].za,ncl[i].lev,reassign);
//    ncl[i].ndisc = getDiscreteLevels(&ncl[i].za,ncl[i].lev,zaid_index,
//		        elev_s_t,p_ng,nf_t,pe_t,ic_t,reassign2);
    statFixDiscreteLevels(&ncl[i]);
    statSetupEnergyBin(&ncl[i]);
    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
//    statClearPopulation(&ncl[i]);
  }

  /*** GDR parameters */
  GDR gdr[MAX_GDR];
  statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,0.0);
  
  /*** special setting for evaporation calc. */
  if(EVAPORATION_SPECTRUM){
    initJ = 0.0;
    initP = 1;
    for(int i=0 ; i<=nemit ; i++){
      ncl[i].ndisc = 1;
      ncl[i].ncont = ncl[i].ntotal;
    }
  }
  
  /*** clear spectrum array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    for(int i=0 ; i<SPECTRA_OUTPUT ; i++) spc[i][k] = 0.0;
  }

  /*** start Monte Carlo */
  if(nsim == 0) nsim = MAX_SIMULATION;

  /*** decay start at a given discrete level */
  if(ncl[0].ncont == 0){
    /*** find nearest discrete level for given Ex */
    int nstart = statStoreLevelExcite(&ncl[0]);

    for(unsigned long t=0 ; t<nsim ; t++){
      np = 0;
      mcGammaCascade(0,nstart,spc);
      if(ctl.print_history) mcWriteHistory();
      recordEmittedParticles();
    }
  }
  /*** start with the continuum bin */
  else{
    /*** set up initial population */
    if(ctl.init_pop == TRANSMISSION){
      CrossSection cx;
      double mu = ncl[1].mass * pdt[neutron].mass / (ncl[1].mass + pdt[neutron].mass);
      tin.ecms = ncl[0].excitation[0] -  ncl[0].cdt[1].binding_energy;
      if(tin.ecms<0.0) cgmTerminateCode("closed neutron channel");
      tin.lmax = omCalc(tin.ecms, &pdt[neutron], &ncl[1].za, mu, tin.tran, &cx);
    }
    statSetupInitialPopulationSpec(initJ,initP,0,spinf,&tin);
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
      for(int j=0 ; j<MAX_J ; j++) {
        pin[k][j] = ncl[0].pop[k][j];			
      }
    }

    MCPoint minit;
    double  emax[MAX_COMPOUND];
    for(int i=0 ; i<=nemit ; i++) emax[i] = ncl[i].max_energy;

    for(unsigned long t=0 ; t<nsim ; t++){
      /*** choose initial state */
      for(int i=0 ; i<=nemit ; i++) ncl[i].max_energy = emax[i];
      
      if(ctl.init_pop == SINGLE) minit = mcFixedStart(emax[0],pin);
      else                       minit = mcSampleStart(emax[0],pin);
      			
      /*** set CN kinetic energy if provided */
      if(targE > 0.0) cnKE = targE;
      
      /*** perform simulation */
      np = 0;
      mcHistoryLoop(minit,pdt,spc);
      if(ctl.print_history) mcWriteHistory();
      recordEmittedParticles();
    }
  }

  /*** restore spectra array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    spc[0][k] = spc[2][k];
    spc[1][k] = spc[3][k];
  }
  mcNormalizeSpectrum(nsim,spc);

  for(int j=0 ; j<MAX_ENERGY_BIN ; j++) delete [] pin[j];
  delete [] pin;
  delete [] tin.tran;
}


/**********************************************************/
/*      Monte Carlo Neutron and Gamma Emission            */
/**********************************************************/
void mcHistoryLoop(MCPoint minit, Pdata *pdt, double **spc)
{
//---------------------------------------
//      Monte Carlo Cascade Calculation
  MCPoint mcurr,mnext;

  mcurr = minit;
  
  while(true){
    /*** prepare energy bins and level densities */
    int c0 = mcurr.getN();
    int c1 = ncl[mcurr.getN()].cdt[neutron].next;

    /*** restore initial population in the top bin,
         and reconstruct bins and level densities */
    mcRestorePopulation(mcurr,&ncl[c0],&ncl[c1]);
		
    /*** calculate Hauser-Feshbach decay */
    if(EVAPORATION_SPECTRUM) specEvaporationBin(c0,pdt,spc);
    else                     spectraBinary(c0,pdt,spc);
        
    /*** move current position to the next state */
    //mnext = mcFindNext(c0,c1);
    mnext = mcFindNext(c0,c1,mcurr.getJ(),mcurr.getP());
		
    if(PERTURB_EXCITATON_ENERGY && !mnext.ifDiscrete()){
      int c = mnext.getN();
      int k = mnext.getK();
      double emax = (k == 0) ? mnext.getE() : (ncl[c].excitation[k-1] +  ncl[c].excitation[k])*0.5;
      double emin = (ncl[c].excitation[k+1] +  ncl[c].excitation[k])*0.5;
      if(emin < ncl[c].lev[ncl[c].ndisc-1].energy) emin = ncl[c].lev[ncl[c].ndisc-1].energy;

			/*** perturbation to the final state energy within its energy bin */
			double ex = emin + (emax-emin)*rng_cgm();
      mnext.setE(ex);
    }

    /*** terminate gamma-decay, if too low continuum energy */

    // ***** FROM IONEL ON SEP. 5, 2012
////    if( (CONTINUUM_LOWER_CUT > 0.0) && (-1 == c1) ){
		if( (CONTINUUM_LOWER_CUT > 0.0) && (c0 == mnext.getN()) ){
			double elcut = ncl[c0].lev[ncl[c0].ndisc-1].energy + CONTINUUM_LOWER_CUT;
      if( (mnext.getE() < elcut) && (!mnext.ifDiscrete()) ){
				int index0=-1 ;
        int ndisc = ncl[c0].ndisc ;
        double xpr = -1.;
        for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
          if( ncl[c0].lev[kk2].parity * mcurr.getP() > 0 ) continue ; //E1 selection rules only
          if( fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() ) < 2. && 1. <= ncl[c0].lev[kk2].spin + mcurr.getJ() ){
            double xpr1 = rng_cgm() ;
            if( xpr1 > xpr ){
              xpr = xpr1 ;
              index0 = kk2;
            }
          }
        }
				
        /* if such a transition does not exist, force transition to the closest-j state */
        if(index0 == -1 ){
          double delJ=100.;
          for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
            if( fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() ) < delJ ){
              delJ = fabs( ncl[c0].lev[kk2].spin - mcurr.getJ() );
              index0 = kk2 ;
            }
          }
          //	  cout << " finished, index0 = " << index0 << endl;
        }
        mnext.set(c0, (int)ncl[c0].lev[index0].parity,(int)ncl[c0].lev[index0].spin,index0,ncl[c0].lev[index0].energy);
        mnext.setDiscrete();
      }
    }

    //#endif
    mcStoreHistory(false,mcurr,mnext,spc);
    mcurr = mnext;

    /*** if discrete levels reached */
    if( mnext.ifDiscrete() ) break;

  }

  /*** discrete level case */
	int final_state = 0;
	if(mnext.getK() > 0) final_state = mcGammaCascade(mnext.getN(),mnext.getK(),spc);

	/*** insert final state */
	int c = mcurr.getN();
	mcurr.set(c,ncl[c].lev[final_state].parity,(int)ncl[c].lev[final_state].spin,final_state,ncl[c].lev[final_state].energy);
	mcurr.setDiscrete();
	hs[np++].set(gammaray, 0.0, mcurr);
	
}


/**********************************************************/
/*      Check Integrant Becomes Larger than Rand          */
/**********************************************************/
MCPoint mcFindNext(int c0, int c1, int jspin, int parity)
{
  MCPoint m;
  double  r = rng_cgm();
  double  s = 0.0;
	
	/* gamma transition to continuum */
  /*** skip k=0, because there are initial populations */
  for(int k=1 ; k<ncl[c0].ncont ; k++){
    for(int j=0 ; j<=ncl[c0].jmax ; j++){
      s += ncl[c0].pop[k][j].even;
      if(s>=r){
        m.set(c0, 1,j,k,ncl[c0].excitation[k]);
        return(m);
      }
      s += ncl[c0].pop[k][j].odd;
      if(s>=r){
        m.set(c0,-1,j,k,ncl[c0].excitation[k]);
        return(m);
      }
    }
  }
	
	/* gamma transition to discrete */
  for(int k=0 ; k<ncl[c0].ndisc ; k++){
    s += ncl[c0].lpop[k];
    if(s>=r){
      m.set(c0,ncl[c0].lev[k].parity,(int)ncl[c0].lev[k].spin,k,ncl[c0].lev[k].energy);
      m.setDiscrete();
      return(m);
    }
  }
	
	/* neutron transition to continuum */
  if(ncl[c0].cdt[neutron].status){
    for(int k=0 ; k<ncl[c1].ncont ; k++){
      for(int j=0 ; j<=ncl[c1].jmax ; j++) {
        s += ncl[c1].pop[k][j].even;
        if(s>=r){
          m.set(c1, 1,j,k,ncl[c1].excitation[k]);
          return(m);
        }
        s += ncl[c1].pop[k][j].odd;
        if(s>=r){
          m.set(c1,-1,j,k,ncl[c1].excitation[k]);
          return(m);
        }
      }
    }
		
		/* neutron transition to discrete */
    for(int k=0 ; k<ncl[c1].ndisc ; k++){
      s += ncl[c1].lpop[k];
      if(s>=r){
        m.set(c1,ncl[c1].lev[k].parity,(int)ncl[c1].lev[k].spin,k,ncl[c1].lev[k].energy);
        m.setDiscrete();
        return(m);
      }
    }
  }

  /*** if no transition, force decay to the ground state */
//  m.set(c0,0,0,0,0.0);
//  //IS: Ionel :: modified to force the transition into a state allowed by multiplicity 1 transition
  int index0=-1 ;
  int ndisc = ncl[c0].ndisc ;
  double xpr = -1.;
  for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
    if( ncl[c0].lev[kk2].parity * parity > 0 ) continue ; //E1 selection rules only
    if( fabs( ncl[c0].lev[kk2].spin - jspin ) < 2. && 1. <= ncl[c0].lev[kk2].spin + jspin ){
      double xpr1 = rng_cgm() ;
      if( xpr1 > xpr ){
        xpr = xpr1 ;
        index0 = kk2;
      }
    }
  }
  /* if such a transition does not exist, force transition to the closest-j state */
  if(index0 == -1 ){
    double delJ=100.;
    for( int kk2 = 0 ; kk2 <  ndisc ; kk2++ ){
      if( fabs( ncl[c0].lev[kk2].spin - jspin ) < delJ ){
        delJ = fabs( ncl[c0].lev[kk2].spin - jspin );
        index0 = kk2 ;
      }
    }
  }
	
  m.set(c0, (int)ncl[c0].lev[index0].parity,(int)ncl[c0].lev[index0].spin,index0,ncl[c0].lev[index0].energy);
  m.setDiscrete();
  
  return(m);
}


/**********************************************************/
/*      Monte Carlo Cascading From Each Level             */
/**********************************************************/
int mcGammaCascade(int c0, int i0, double **spc)
{
  MCPoint mcurr,mnext;

  int m = 0;
  double tstamp=0.;
  int igam=0;

  while(1){

    mcurr.set(c0,ncl[c0].lev[i0].parity,(int)ncl[c0].lev[i0].spin,i0,ncl[c0].lev[i0].energy);
    mcurr.setDiscrete();
    if(i0==0) return 0;

    /*
    if (ncl[c0].lev[i0].halflife>0.0 && EXPERIMENTAL_TIME_WINDOW>0.0) {
      if (rng_cgm()>1.0-exp(-EXPERIMENTAL_TIME_WINDOW/ncl[c0].lev[i0].halflife*LOG_2)) return(i0);
    }

    if (ncl[c0].lev[i0].halflife>0.0) {
      double trnd = 10.* ncl[c0].lev[i0].halflife * rng_cgm() ; // sample a time of up to 10 times the lifetime
      double tau = ncl[c0].lev[i0].halflife / LOG_2;
      while( rng_cgm() > exp(-trnd/tau) ){
        trnd = 10.* ncl[c0].lev[i0].halflife * rng_cgm() ;
      }
      tstamp += trnd ;
    }

    */

    // IS: changed the sampling of the emission time

    if (ncl[c0].lev[i0].halflife>0.0) {
      tstamp += -log(rng_cgm())*ncl[c0].lev[i0].halflife / LOG_2;
    }
    if(time_coincidence_window>0.0 && tstamp>time_coincidence_window)
      break;
    

    int zaid = ncl[c0].za.getZ()*1000+ncl[c0].za.getA();
    double r = rng_cgm();
    double s = 0.0;

    for(int j=0 ; j<ncl[c0].lev[i0].ngamma ; j++){
      int i1 = ncl[c0].lev[i0].fstate[j];
      s += ncl[c0].lev[i0].branch[j];
      if(s >= r){
        /*** if Internal Conversion happens, store no gamma information */
        bool ic = false;
        if( INCLUDE_INTERNAL_CONVERSION && (rng_cgm() > ncl[c0].lev[i0].gratio[j]) )
	  ic = true;
	gammaTimeStamp[np]=tstamp;
        /*** point next state */
        mnext.set(c0,ncl[c0].lev[i1].parity,(int)ncl[c0].lev[i1].spin,i1,ncl[c0].lev[i1].energy);
        mnext.setDiscrete();
        /*** save history */
        mcStoreHistory(ic,mcurr,mnext,spc);
        m++;
        i0 = i1;
        mcurr = mnext;
        break;
      }
    }

    if(i0 == 0 ) break;
  }
  return(i0);
}


/**********************************************************/
/*      Sample JP from Initial Distribution               */
/**********************************************************/
MCPoint mcSampleStart(double e0, Parity **p)
{
  int     n0 = 0, k0 = 0;
  MCPoint m(n0,0,0,k0,e0);

  double  x = 0.0;
  double  r = rng_cgm();
  for(int j=0 ; j<MAX_J ; j++){
    x += p[k0][j].even;
    if(x >= r){ m.set(n0, 1,j,k0,e0); return(m); }
    x += p[k0][j].odd ;
    if(x >= r){ m.set(n0,-1,j,k0,e0); return(m); }
  }

  m.set(n0,1,MAX_J-1,k0,e0);
  return(m);
}


/**********************************************************/
/*      Search for Starting JP State                      */
/**********************************************************/
MCPoint mcFixedStart(double e0, Parity **p)
{
  int     n0 = 0, k0 = 0;
  MCPoint m(n0,0,0,k0,e0);

  /*** starting point, if fixed */
  bool found = false;
  for(int j=0 ; j<MAX_J ; j++){
    if(     p[k0][j].even > 0.0){
      m.set(n0, 1,j,k0,e0); found = true; break;
    }
    else if(p[k0][j].odd  > 0.0){
      m.set(n0,-1,j,k0,e0); found = true; break;
    }
  }
  if(!found) cgmTerminateCode("initial population zero");

  return(m);
}


/**********************************************************/
/*      Store Initial Population in the Top Bin           */
/**********************************************************/
void mcRestorePopulation(MCPoint m, Nucleus *n0, Nucleus *n1)
{
  n0->max_energy = m.getE();
  statSetupEnergyBin(n0);
  statSetupLevelDensity(n0,&n0->ldp);
//  statClearPopulation(n0);

  if(n0->cdt[neutron].status){
    n1->max_energy = n0->max_energy - n0->cdt[neutron].binding_energy;
    statSetupEnergyBin(n1);
    statSetupLevelDensity(n1,&n1->ldp);
//    statClearPopulation(n1);
  }

  if(m.getP() > 0) n0->pop[0][m.getJ()].even = 1.0;
  else             n0->pop[0][m.getJ()].odd  = 1.0;
}


/**********************************************************/
/*      Store Decay Gamma-Ray Energies                    */
/**********************************************************/
void mcStoreHistory(bool ic, MCPoint m0, MCPoint m1, double **spc)
{
  double ecms = 0.0;
  int k = 0;
  
  if(m0.getN() == m1.getN()){
    ecms = m0.getE() - m1.getE();
    if(ic) hs[np++].set(unknown , ecms, m0);
    else   hs[np++].set(gammaray, ecms, m0);
  }else{
    /*** if neutron, subtract separation energy */
    ecms = m0.getE() - m1.getE() - ncl[m0.getN()].cdt[neutron].binding_energy;
    hs[np++].set(neutron, ecms, m0);
  }
  if(np == MAX_CASCADE) cgmTerminateCode("too many particle decays");
  
#ifdef HAVE_PRIVATE_ENERGY_GRID
  k = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,ecms,ncl[0].binwidth);
#else
  k = specFindEnergyBin(ecms,ncl[0].de);
#endif

  if(k > 0){
    if(m0.getN() == m1.getN()) spc[2][k] += 1.0;
    else                       spc[3][k] += 1.0;
  }
}





/**********************************************************/
/*      Write Decay History                               */
/**********************************************************/
void mcWriteHistory()
{
  double  etot = 0.0;
  char    name[9] = {'g','n',' ',' ',' ',' ',' ',' ','e'};
  
  /*** convert neutron energy into lab-energy */
  if( cnKE > 0.0 ){
    for(int i=0 ; i<np-1 ; i++){
      if(hs[i].getIndex() == 1){
        MCPoint mp = hs[i].getMCPoint();
        double  ac = ncl[ mp.getN() ].za.getA();
        double  ec = hs[i].getEnergy();
        double  el = mcCmsToLabEnergy(ec,ac);
        hs[i].set(neutron,el,mp);
      }
    }
  }
  
  /*** the last element of hs[] contains ground state spin and parity only */
  if(EVENT_OUTPUT_FORMAT == 1){
    int j = np-1;
    for(int i=0 ; i<np-1 ; i++) if(hs[i].getIndex() == unknown) j--;

    cout << setw(5) << j << endl;
    cout << setprecision(5) << setiosflags(ios::scientific);

    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      if(hs[i].getIndex() == unknown) continue;
      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    cout << endl;
  }

  else if(EVENT_OUTPUT_FORMAT == 2){
    cout << setprecision(5) << setiosflags(ios::scientific);

    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      cout << setw(4)  << name[hs[i].getIndex()];
      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    cout << endl;

    /*** print out parent state */
    cout << "     ";
    for(int i=0 ; i<np ; i++){
      MCPoint m = hs[i].getMCPoint();
      cout << setw(3)  << m.getJ();
      cout << setw(1)  << ((m.getP() == -1) ? '-' : '+');
      cout << setw(13) << m.getE();
    }
    cout << endl;
  }

  else{
    //    cout << setprecision(5) << setiosflags(ios::scientific);
    //    cout << setw(5) << np-1;
    
    /*** print out particle energies */
    for(int i=0 ; i<np-1 ; i++){
      //      cout << setw(4)  << name[hs[i].getIndex()];
      //      cout << setw(13) << hs[i].getEnergy();
      etot += hs[i].getEnergy();
    }
    //  cout << " : " << etot;
    //    cout << endl;
  }
	
}


/**********************************************************/
/*      Normalize MC Spectrum                             */
/**********************************************************/
void mcNormalizeSpectrum(unsigned long n, double **spc)
{
  for(int i=0 ; i<2 ; i++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[i][k] = spc[i][k]/n;
  }
}


/**********************************************************/
/*      Determine Lab Energy of Emitted Neutron           */
/**********************************************************/
double mcCmsToLabEnergy(double ecms, double mass)
{
  /*** sample emission angle */
  double costh = cos(PI*rng_cgm());
  double elab  = ecms + cnKE/mass + 2.0*sqrt(ecms*cnKE/mass)*costh;
  
  /*** set new recoil energy */
  double cosph = sqrt(1.0 - ecms/elab * (1.0 - costh*costh));
  cnKE = ( elab + mass*cnKE - 2.0*sqrt(elab*cnKE*mass)*cosph )/(mass-1.0);
  //cnKE = ecms/(mass-1.0)+(mass-1.0)/mass*cnKE - 2.*sqrt(ecms*cnKE/mass)*costh;
  //cout << ecms <<"  "<< elab <<"  :" << acos(costh) / PI * 180<<endl;
  
  return(elab);
}


/**********************************************************/
/*      For Monitoring Only                               */
/**********************************************************/
#ifdef DEBUG_POPCHECK
#include <iomanip>
void specPopCheck(Nucleus *n)
{
  double s1 = 0.0, s2 = 0.0;

  cout << "#  ";
  cout << setw(3) << setfill('0') << n->za.getZ() << '-';
  cout << setw(3) << setfill('0') << n->za.getA() << endl;
  cout << setfill(' ');

  for(int k=0 ; k<n->ncont ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(7) << setprecision(3) <<  n->excitation[k];
    cout.setf(ios::scientific, ios::floatfield);

    double s3 = 0.0;
    for(int j=0 ; j<19 ; j++) {
      s1 += n->pop[k][j].even+n->pop[k][j].odd;
      s3 += n->pop[k][j].even+n->pop[k][j].odd;
      cout << setprecision(2) << setw(9) << n->pop[k][j].even+n->pop[k][j].odd;
    }
    cout << setprecision(2) << setw(9) << s3 << endl;
  }

  s2 = 0.0;
  for(int k=0;k<n->ndisc;k++){
    s2 += n->lpop[k];
    printf("%5d %10.6f %10.3e\n",k,n->lev[k].energy,n->lpop[k]);
  }

  cout << "# SUMcont " << s1 << endl;
  cout << "# SUMdisc " << s2 << endl;
  cout << "# SUMall  " << s1+s2 << endl;
  cout << endl;
  cout << endl;
}
#endif


/********************************************************************************
 CGMF routine to record emitted neutrons and photons from fission event.
 *******************************************************************************/
void recordEmittedParticles (void) {

	int index;
	MCPoint mc_curr, mc_next;
	
  double neutronEnergies [MAX_NUMBER_PARTICLES];
  double gammaEnergies [MAX_NUMBER_PARTICLES];
	int gammaTypes [MAX_NUMBER_PARTICLES];
	double icEnergies [MAX_NUMBER_PARTICLES];
	
  int neutronMultiplicity=0;
  int gammaMultiplicity=0;
	int icMultiplicity=0;
	
	std::fill_n (neutronEnergies, MAX_NUMBER_PARTICLES, 0.0);
	std::fill_n (gammaEnergies, MAX_NUMBER_PARTICLES, 0.0);
		
  MCPoint mp0 = hs[0].getMCPoint();
  int A = ncl[mp0.getN()].za.getA();
	int Z = ncl[mp0.getN()].za.getZ();
	
	float oddSpin=0.0; if (A%2) { oddSpin=0.5; }
	
	float  Ji = mp0.getJ() + oddSpin; // retrieve initial spin index
	int    Pi = mp0.getP(); // retrieve initial parity index
	double Ui = mp0.getE(); // retrieve initial excitation energy

	for (int i=0; i<np-1; i++) { // loop over emitted particles
		index = hs[i].getIndex();
		mc_curr = hs[i].getMCPoint();
		mc_next = hs[i+1].getMCPoint();
		if (index==1) { // neutrons
			neutronEnergies[neutronMultiplicity++] = hs[i].getEnergy();
		} else if (index==0) { // gammas (no IC)
			gammaTypes[gammaMultiplicity] = 1; // continuum-to-continuum
			if (mc_curr.ifDiscrete()) {
				gammaTypes[gammaMultiplicity] = 3; // discrete-to-discrete
			} else if (mc_next.ifDiscrete()) {
				gammaTypes[gammaMultiplicity] = 2; // continuum-to-discrete
			}
			gammaEnergies[gammaMultiplicity++] = hs[i].getEnergy();
		} else { // Internal Conversion
			icEnergies[icMultiplicity++] = hs[i].getEnergy();
		}
	}
  
}

/********************************************************************************
 CGMF routine to record emitted neutrons and photons from fission event.
 *******************************************************************************/
void recordEmittedParticlesFission (fissionEventType *fe) {
	
	int index;
	MCPoint mc_curr, mc_next;
	
	double cmNeutronEnergies  [MAX_NUMBER_NEUTRONS];
	double neutronEnergies	  [MAX_NUMBER_NEUTRONS];
	double cmGammaEnergies    [MAX_NUMBER_GAMMAS];
	double gammaEnergies      [MAX_NUMBER_GAMMAS];
	double timeStamp          [MAX_NUMBER_GAMMAS];
	double icEnergies		      [MAX_NUMBER_PARTICLES];
	
	int neutronMultiplicity=0;
	int gammaMultiplicity=0;
	int icMultiplicity=0;
	
	std::fill_n (neutronEnergies, MAX_NUMBER_NEUTRONS, 0.0);
	std::fill_n (gammaEnergies, MAX_NUMBER_GAMMAS, 0.0);
	std::fill_n (timeStamp, MAX_NUMBER_GAMMAS, 0.0);
	
  double fragmentMomentum [3] = {};
  for (int i=0; i<3; i++) fragmentMomentum[i]=fe->preFragmentMomentum[i];

  MCPoint mp0 = hs[0].getMCPoint();
  int A = ncl[mp0.getN()].za.getA();
  int Z = ncl[mp0.getN()].za.getZ();
  
  float oddSpin=0.0; if (A%2) { oddSpin=0.5; }
  
  for (int i=0; i<np-1; i++) { // loop over emitted particles
    index = hs[i].getIndex();
    mc_curr = hs[i].getMCPoint();
    mc_next = hs[i+1].getMCPoint();
    if (index==1) { // neutrons
      cmNeutronEnergies[neutronMultiplicity++] = hs[i].getEnergy();
    } else if (index==0) { // gammas (no IC)
      if (mc_curr.ifDiscrete()) {
	timeStamp[gammaMultiplicity] = gammaTimeStamp[i]; // discrete-to-discrete
      }
      cmGammaEnergies[gammaMultiplicity++] = hs[i].getEnergy();
    } else { // Internal Conversion
      icEnergies[icMultiplicity++] = hs[i].getEnergy();
    }
  }
	
	fe->nu  = neutronMultiplicity;
	fe->nug = gammaMultiplicity;
	
	if (neutronMultiplicity>MAX_NUMBER_NEUTRONS) cgmTerminateCode ("[specmc] Too large neutron multiplicity!");
	if (gammaMultiplicity>MAX_NUMBER_GAMMAS) cgmTerminateCode ("[specmc] Too large gamma multiplicity");
	
	double Vn[MAX_NUMBER_NEUTRONS][3];
	double Vg[MAX_NUMBER_GAMMAS][3];
	double cmVn[MAX_NUMBER_NEUTRONS][3];
	double cmVg[MAX_NUMBER_GAMMAS][3];
	double vx,vy,vz,v;
	
	// from kinematics.cpp/.h
  for (int i=0; i<3; i++) fe->postFragmentMomentum[i] = fe->preFragmentMomentum[i];
	boost (fragmentMomentum, Z, A,
				 neutronMultiplicity, cmNeutronEnergies, neutronEnergies, cmVn, Vn,
				 gammaMultiplicity, cmGammaEnergies, gammaEnergies, cmVg, Vg);
	for (int i=0; i<3; i++) fe->postFragmentMomentum[i] = fragmentMomentum[i];

  // post-neutron emission fragment kinetic energy
  double mass = A*amuMeV+getMassExcess(Z*1000+A);
  for (int i=0; i<3; i++) fe->KEpost += ( fe->postFragmentMomentum[i]*fe->postFragmentMomentum[i] );
  fe->KEpost /= (2.0*mass);

	for (int i=0; i<fe->nu; i++) {
		fe->neutronEnergies[i]   = neutronEnergies[i];
		vx = Vn[i][0];
		vy = Vn[i][1];
		vz = Vn[i][2];
		v = sqrt(vx*vx+vy*vy+vz*vz);
		fe->neutronDircosu[i]    = vx/v;
		fe->neutronDircosv[i]    = vy/v;
		fe->neutronDircosw[i]    = vz/v;
		
		fe->cmNeutronEnergies[i]   = cmNeutronEnergies[i];
		vx = cmVn[i][0];
		vy = cmVn[i][1];
		vz = cmVn[i][2];
		v = sqrt(vx*vx+vy*vy+vz*vz);
		fe->cmNeutronDircosu[i]    = vx/v;
		fe->cmNeutronDircosv[i]    = vy/v;
		fe->cmNeutronDircosw[i]    = vz/v;
		
	}
	
	int igc=0; // count the continuum-to-continuum and continuum-to-discrete transitions
	for (int i=0; i<fe->nug; i++) {
	  fe->photonAges[i] = timeStamp[i];
		fe->photonEnergies[i]   = gammaEnergies[i];
		vx = Vg[i][0];
		vy = Vg[i][1];
		vz = Vg[i][2];
		v = sqrt(vx*vx+vy*vy+vz*vz);
		//		fe->photonVelocities[i] = v;
		fe->photonDircosu[i]    = vx/v;
		fe->photonDircosv[i]    = vy/v;
		fe->photonDircosw[i]    = vz/v;
		
		fe->cmPhotonEnergies[i] = cmGammaEnergies[i];
		vx = cmVg[i][0];
		vy = cmVg[i][1];
		vz = cmVg[i][2];
		v = sqrt(vx*vx+vy*vy+vz*vz);
		//		fe->photonVelocities[i] = v;
		fe->cmPhotonDircosu[i]    = vx/v;
		fe->cmPhotonDircosv[i]    = vy/v;
		fe->cmPhotonDircosw[i]    = vz/v;
	}
	
	return;
	
}

