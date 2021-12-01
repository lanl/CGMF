/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf_help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file specmain.cpp

  \brief Main particle emission spectrum calculation

*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <functional>

using namespace std;

#include "cgm.h"
#include "global.h"
#include "ripl2levels.h"
#include "terminate.h"
#include "config.h"
#include "rngcgm.h"

static void    specCompoundDecay   (int, int, Transmission *, Transmission *, double **, double **);
static double  specTransmissionSum (Statcalcmode, bool *, int, int, int, double **,
                                    Transmission *, Transmission *, double, Nucleus *, double **);
static void    specPopulationBypass(int k0, double dp, Nucleus *n, double *spc, int jspin , int parity );

#undef DEBUG_POPCHECK
#ifdef DEBUG_POPCHECK
static void    specPopCheck        (Nucleus *);
#endif

/**********************************************************/
/*      Main Spectra Calculation                          */
/**********************************************************/
void specMain(double initJ, int initP, double spinf, double **spc, double beta2=0.0)
{
  Transmission  tin;
  Pdata pdt[MAX_CHANNEL];
    
#ifdef  HAVE_PRIVATE_ENERGY_GRID
  if(BINARY_REACTION_ONLY == false){
    cgmTerminateCode("Private energy grid should be for binary reaction");
  }
#endif
  
  try{
    tin.tran = new double [3*MAX_J];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in specMain");
  }
  
  //---------------------------------------
  //      Parameter Setting
  
  /*** how many neutrons can be emitted */
  int nemit = 0;
  if(INCLUDE_NEUTRON_EMISSION){
    nemit = statTotalNeutrons(ncl[0].max_energy,&ncl[0].za);
  }
  if(BINARY_REACTION_ONLY){
    if(nemit > 1) nemit = 1;
  }
  
  
  /*** initialize system parameters */
  statSetupInitSystem(nemit,pdt);
	getRiplDiscreteLevels (ncl, nemit);
  
  /*** setup discrete levels, energy bins and level densities */
  for(int i=0 ; i<=nemit ; i++){
    statFixDiscreteLevels(&ncl[i]);
    statSetupEnergyBin(&ncl[i]);
    statSetupLevelDensityParameter(&ncl[i],&ncl[i].ldp);
    statClearPopulation(&ncl[i]);
  }
  
  /*** GDR parameters */
  GDR gdr[MAX_GDR];
  statSetupGdrParameter(&ncl[0],ncl[0].ldp.a,gdr,beta2);
  
  /*** store initial population in the top bin */
  int nstart = 1;
  if(ncl[0].ncont > 0){
    
    if(ctl.init_pop == TRANSMISSION){
      CrossSection cx;
      double mu = ncl[1].mass * pdt[neutron].mass / (ncl[1].mass + pdt[neutron].mass);
      tin.ecms = ncl[0].excitation[0] -  ncl[0].cdt[1].binding_energy;
      if(tin.ecms<0.0) cgmTerminateCode("closed neutron channel");
      tin.lmax = omCalc(tin.ecms, &pdt[neutron], &ncl[1].za, mu, tin.tran, &cx);
    }
    
    statSetupInitialPopulationSpec(initJ,initP,0,spinf,&tin);
    nstart = ncl[0].ndisc;
  }
  /*** find nearest discrete level for given Ex */
  else nstart = statStoreLevelExcite(&ncl[0]);
  
  /*** Hauser-Feshbach calculation */
  if(BINARY_REACTION_ONLY){
    if(ncl[0].ncont > 0) spectraBinary(0,pdt,spc);
    else                 specGammaDiscrete(nstart,spc[0],&ncl[0]);
  }else if(EVAPORATION_SPECTRUM){
    if(ncl[0].ncont > 0) specEvaporation(nemit,pdt,spc);
  }else{
    if(ncl[0].ncont > 0) spectra(nemit+1,pdt,spc);
    else                 specGammaCascade(nstart,spc[0],&ncl[0]);
  }
  
  delete [] tin.tran;
  
  return;
}


/**********************************************************/
/*      Gamma and Neutron Emission Spectrum Calculation   */
/**********************************************************/
void spectra(int nc, Pdata *pdt, double **spc)
{
  Transmission  *tc = NULL;
  Transmission  *td = NULL;
  double        *tg[MAX_MULTIPOL], *gn[2];
  
  //---------------------------------------
  //      Allocate Transmission Data
  
  try{
    tc = new Transmission [MAX_ENERGY_BIN];
    td = new Transmission [MAX_LEVELS    ];
    
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
    for(int k=0 ; k<MAX_LEVELS     ; k++) td[k].tran = new double [3*MAX_J];
    for(int i=0 ; i<MAX_MULTIPOL ; i++)   tg[i]      = new double [MAX_ENERGY_BIN];
    for(int i=0 ; i<2 ; i++)              gn[i]      = new double [MAX_ENERGY_BIN];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in spectra");
  }
  
  /*** clear total spectrum array */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[0][k] = spc[1][k] = spc[2][k] = spc[3][k] = 0.0;
  
  //---------------------------------------
  //      Cascade Calculation
  
  /*** for all CN */
  for(int c0=0; c0<nc ; c0++){
    
    /*** clear neutron and gamma-ray spectrum array */
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) gn[0][k] = gn[1][k] = 0.0;
    
    /*** calculate transmission coefficients for neutron channel */
    if(ncl[c0].cdt[neutron].status) statStoreContinuumTransmission(c0,0,&pdt[neutron],tc);
    
    /*** loop over parent excitation */
    for(int k0=0 ; k0<ncl[c0].ncont ; k0++){
      
      /*** store gamma-ray transmission data into array */
      statStoreGammaTransmission(k0,tg,&ncl[c0]);
      if(ncl[c0].cdt[neutron].status) statStoreDiscreteTransmission(c0,k0,&pdt[neutron],td);
      
      /*** CN decay */
      specCompoundDecay(c0,k0,tc,td,tg,gn);
      
#ifdef DEBUG_POPCHECK
      if(k0 == 0) specPopCheck(&ncl[c0]);
#endif
    }
    
    /*** gamma-ray cascading, ground state population */
    specGammaCascade(ncl[c0].ndisc,gn[0],&ncl[c0]);
    
    /*** cumulative neutron and gamma-ray spectra */
    for(int k=0 ; k<ncl[c0].ntotal ; k++){
      for(int i=0 ; i<2 ; i++) spc[i][k] += gn[i][k];
    }
    
    if(ctl.print_indivspec) cgmPrintSpectra(false,ncl[0].de,gn);
  }
  
  
  //---------------------------------------
  //      Free Allocated Memory
  
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[k].tran;
  for(int k=0 ; k<MAX_LEVELS     ; k++) delete [] td[k].tran;
  for(int i=0 ; i<MAX_MULTIPOL ; i++)   delete [] tg[i];
  for(int i=0 ; i<2 ; i++)              delete [] gn[i];
  
  delete [] tc;
  delete [] td;
}


/**********************************************************/
/*      Binary Reaction Calculation                       */
/**********************************************************/
void spectraBinary(int c0, Pdata *pdt, double **spc)
{
  Transmission  *tc = NULL;
  Transmission  *td = NULL;
  double        *tg[MAX_MULTIPOL];
  
  //---------------------------------------
  //      Allocate Transmission Data
  
  try{
    tc = new Transmission [MAX_ENERGY_BIN];
    td = new Transmission [MAX_LEVELS    ];
    
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[k].tran = new double [3*MAX_J];
    for(int k=0 ; k<MAX_LEVELS     ; k++) td[k].tran = new double [3*MAX_J];
    for(int i=0 ; i<MAX_MULTIPOL ; i++)   tg[i]      = new double [MAX_ENERGY_BIN];
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in spectra");
  }
  
  /*** clear total spectrum and energy arrays */
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[0][k] = spc[1][k] = 0.0;
  
  
  //---------------------------------------
  //      Binary Reaction Calculation
  
  if(ncl[c0].cdt[neutron].status){
    statStoreContinuumTransmission(c0,0,&pdt[neutron],tc);
    statStoreDiscreteTransmission(c0,0,&pdt[neutron],td);
  }
  statStoreGammaTransmission(0,tg,&ncl[c0]);
  
  specCompoundDecay(c0,0,tc,td,tg,spc);
  
  
  //---------------------------------------
  //      Free Allocated Memory
  
  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[k].tran;
  for(int k=0 ; k<MAX_LEVELS     ; k++) delete [] td[k].tran;
  for(int i=0 ; i<MAX_MULTIPOL ; i++)   delete [] tg[i];
  
  delete [] tc;
  delete [] td;
  
#ifdef DEBUG_POPCHECK
  specPopCheck(&ncl[c0]);
  if(ncl[c0].cdt[neutron].status){
    specPopCheck(&ncl[ ncl[c0].cdt[neutron].next ]);
  }
#endif
}


/**********************************************************/
/*      Particle Emission from CNparent to CNdaughter     */
/*      -----                                             */
/*             Usual Hauser-Feshhach calculation          */
/**********************************************************/
void    specCompoundDecay(int c0, int k0, Transmission *tc, Transmission *td,
                          double **tg, double **spc)
{
  bool   tstat[MAX_CHANNEL];
  int    i0 = (int)(2.0*halfint(ncl[c0].lev[0].spin));
  
  /*** Check if channel is closed */
  tstat[0] = true;
  tstat[1] = ncl[c0].cdt[neutron].status;
  if(tstat[1] && (ncl[c0].excitation[k0] < ncl[c0].cdt[neutron].binding_energy)){
    tstat[1] = false;
  }
  
  /*** skip gamma-neutron channel if Ex > Bn + IGNORE_CAPTURE_GAMMA_RAY */
  if(IGNORE_CAPTURE_GAMMA_RAY > 0.0){
    if( (c0 == 0) && tstat[1] ){
      double excut = ncl[c0].cdt[neutron].binding_energy + IGNORE_CAPTURE_GAMMA_RAY;
      if(ncl[c0].excitation[k0] > excut ) tstat[0] = false;
    }
  }
  
  /*** Loop over CN J and Pariry */
  for(int j0=i0 ; j0<MAX_J*2+i0 ; j0+=2){
    int jdx = (j0-i0)/2;
    for(int p0=-1 ; p0<=1 ; p0+=2){
      
      double pop = (p0==1) ? ncl[c0].pop[k0][jdx].even : ncl[c0].pop[k0][jdx].odd;
      if( pop==0.0 ) continue;
      
      double tsum=0.0;
      
      /*** Sum over gamma and particle decay */
      tsum  = specTransmissionSum(sumall,tstat,k0,p0,j0,tg,tc,td,0.0,&ncl[c0],spc);
      
      /*** If no transition, add population to g.s */
      if(tsum==0.0 && pop>0.0){
        specPopulationBypass(k0,pop,&ncl[c0],spc[0], jdx, p0);
        continue;
      }
      
      /*** Population add to each resudual CN */
      double dp = (tsum == 0.0) ? 0.0 : pop/tsum;
      specTransmissionSum(hauser,tstat,k0,p0,j0,tg,tc,td,dp,&ncl[c0],spc);
    }
  }
  
}


/**********************************************************/
/*      Sum All Transmission Coefficients from CNparent   */
/**********************************************************/
double  specTransmissionSum(Statcalcmode mode, bool *tstat, int k0, int pcn, int jcn,
                            double **tg, Transmission *tc, Transmission *td,
                            double dpop, Nucleus *n0, double **spc)
{
  double tsum=0.0;
  
  /*** Sum up gamma-ray transmissions */
  if(tstat[0]) tsum = specTransitionGamma(mode,k0,pcn,jcn,tg,dpop,n0,spc[0]);
  
  /*** Sum up neutron transmissions */
  if(tstat[1]){
    Nucleus *n1 = &ncl[n0->cdt[neutron].next];
    tsum += specTransitionParticle(mode,k0,pcn,jcn,tc,td,dpop,n0,n1,spc[1]);
  }
  
  return(tsum);
}


/**********************************************************/
/*      Avoid Population Trap for No Transition Case      */
/*      -----                                             */
/*             In case J in the continuum is so high,     */
/*             there is no transition to the discrete     */
/*             levels due to low multipolarity.           */
/*             These trapped flux is directy fed to       */
/*             the ground state                           */
/**********************************************************/
/* Ionel, 8/29/2013 */
void    specPopulationBypass(int k0, double dp, Nucleus *n, double *spc,
                             int jspin , int parity )
{
  double eg = n->excitation[k0];
  
  int kk0 = -1 ; double delJ = 1000. ;
  double xpr = -1.;
  for( int kk2 = 0 ; kk2 <  n->ndisc ; kk2++ ){
    if( n->lev[kk2].parity * parity > 0 ) continue ; //E1 selection rules only
    if( fabs( n->lev[kk2].spin - jspin ) < 2. && 1. <= n->lev[kk2].spin + jspin ){
			double xpr1 = rng_cgm();
      if( xpr1 > xpr ){
        xpr = xpr1 ;
        kk0 = kk2;
      }
    }
  }
  /* if such a transition does not exist, force transition to the
   closest-j state */
  if(kk0 == -1 ){
    for( int kk2 = 0 ; kk2 <  n->ndisc ; kk2++ ){
      if( delJ >= fabs( n->lev[kk2].spin - jspin ) ){
        delJ = fabs( n->lev[kk2].spin - jspin ) ;
        kk0=kk2 ;
      }
    }
  }
  eg -= n->lev[kk0].energy ;
  
  
#ifdef HAVE_PRIVATE_ENERGY_GRID
  int kg = specFindCustomGrid(NUMBER_OF_PRIVATE_GRID,eg,n->binwidth);
#else
  int kg = specFindEnergyBin(eg,n->de);
#endif
  
  //  if(kg >= 0) spc[kg] += dp;
  n->lpop[kk0] += dp;
}


/**********************************************************/
/*      Find a Bin Index for Given Excitation             */
/**********************************************************/
int    specFindEnergyBin(double e, double de)
{
  double eps = 1.0e-4;
  
  if(fabs(e)<=eps) return(0);
  if(e<0.0 || e> MAX_ENERGY_BIN*de) return(-1);
  
  if(e<de*0.5) return(0);
  
  /*** Add first energy bin, which has a half width,
   so that zero-energy emission is mapped onto k=0 */
  e += 0.5*de;
  
  /*** To avoid round-off error, a small number is added */
  int k = (int)(e/de + eps );
  
  if( (de*k<= (e+eps)) && (e<de*(k+1)) ) return (k);
  
  cerr << "out of range " << k << " "<<de*k <<  " " << e << " " << de*(k+1) << endl;
  return(-1);
}


/**********************************************************/
/*      Find a Bin Index for Custom Energy Strucure       */
/**********************************************************/
int    specFindCustomGrid(int ng, double e, double *de)
{
  int k=0;
  double emin, emax;
  
  bool inrange = false;
  
  emin = 0.0;
  emax = de[0]*0.5;
  for(k=0 ; k<ng ; k++){
    if( (emin <= e) && (e < emax) ){
      inrange = true;
      break;
    }
    emin = emax;
    emax += de[k+1];
    //  cout << k <<" " << emin <<" "<< emax <<" " << e << endl;
  }
  
  if(inrange) return(k);
  else{
    cerr << "out of range " << e << endl;
    return(-1);
  }
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
