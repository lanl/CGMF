/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file setup.cpp

  \brief Initialize all quantities in nucleus and GDR parameters

*/

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <functional>

// TAW
// added 8/4/2011
#include <time.h>

using namespace std;

#include "cgm.h"
#include "kcksyst.h"
#include "masstable.h"
#include "terminate.h"
#include "config.h"
#include "physics.h"
#include "rngcgm.h"

//#include "global_var.h"
#include "kcksyst.h"
//#include "kcksyst2.h"

#ifdef  HAVE_PRIVATE_ENERGY_GRID
#include GRID_STRUCTURE_FILE
#endif

static inline double binenergy(double *, int);
static inline double normaldist(double, double);

/**********************************************************/
/*     Initialize All Data                                */
/**********************************************************/
void statSetupInitSystem(int nemit, Pdata *pdt)
{
  /*** store particle mass, spin, and identifiers */
  pdt[gammaray].particle.setZA(0,0);
  pdt[gammaray].particleID  = gammaray;
  pdt[gammaray].omp         = 0;
  pdt[gammaray].spin        = 0.0;
  pdt[gammaray].mass        = 0.0;
  pdt[gammaray].mass_excess = 0.0;

  pdt[neutron].particle.setZA(0,1);
  pdt[neutron].particleID   = neutron;
  pdt[neutron].omp          = 6 << 8; // Koning-Delaroche
  pdt[neutron].spin         = 0.5;
  pdt[neutron].mass         = NEUTRONMASS;
  pdt[neutron].mass_excess  = ENEUTRON;

  /*** Z and A number after neutron emission */
  ncl[1].za = ncl[0].za - pdt[neutron].particle;
  if(nemit >1 ){
    for(int i=2 ; i<=nemit ; i++) ncl[i].za = ncl[i-1].za - pdt[neutron].particle;
  }

  for(int i=0 ; i<=nemit ; i++){
    ncl[i].mass_excess = mass_excess(ncl[i].za.getZ(),ncl[i].za.getA());
    ncl[i].mass = ncl[i].za.getA() + ncl[i].mass_excess / AMUNIT;

    for(int j=0 ; j<MAX_CHANNEL ; j++){
      ncl[i].cdt[j].status         = false;
      ncl[i].cdt[j].particleID     = (Particle)j;
      ncl[i].cdt[j].next           =  -1;
      ncl[i].cdt[j].binding_energy = 0.0;
      ncl[i].cdt[j].spin2          = (int)(2.0*pdt[j].spin);
    }
  }

  if(INCLUDE_NEUTRON_EMISSION){
    for(int i=1 ; i<=nemit ; i++){

      double sepa = ncl[i].mass_excess + pdt[neutron].mass_excess - ncl[i-1].mass_excess;
      ncl[i-1].cdt[neutron].binding_energy = sepa;

      if(sepa < ncl[i-1].max_energy){
        ncl[i-1].cdt[neutron].status = true;
        ncl[i-1].cdt[neutron].next   = i;
        ncl[i  ].max_energy          = ncl[i-1].max_energy - sepa;
      }
    }
  }
}


/**********************************************************/
/*     Generate Energy Mesh in the Continuum              */
/*      ----                                              */
/*          When the highest discrete level energy is     */
/*          higher than the maximum excitation (elmax),   */
/*          generate continuum bins.                      */
/**********************************************************/
int statSetupEnergyBin(Nucleus *n)
{
  int num_levels = (n->ndisc == 0) ? 0 : n->ndisc - 1;
  double elmax = n->lev[num_levels].energy;
  int    khigh=MAX_ENERGY_BIN;
  Parity pzero;

  n->ncont = 0;
  n->de    = ENERGY_BIN;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){
    n->excitation[k] = 0.0;

#ifdef  HAVE_PRIVATE_ENERGY_GRID
    if(k == 0)
      n->binwidth[k] = custom_energy_grid[1];
    else if(k < NUMBER_OF_PRIVATE_GRID-1)
      n->binwidth[k] = (custom_energy_grid[k+1] - custom_energy_grid[k-1])*0.5;
    else
      n->binwidth[k] = ENERGY_BIN;
#else
    n->binwidth[k]   = ENERGY_BIN;
#endif
    for(int j=0 ; j<MAX_J ; j++){
      n->pop[k][j]      = pzero;
//      n->density[k][j]  = pzero;
    }
  }
  for(int i=0 ; i<MAX_LEVELS ; i++) n->lpop[i]=0.0;

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++){

    /*** bin energy is the mid-point of each bin,
         Khigh is the bin number of highest excitation */
    if(binenergy(n->binwidth,k) > n->max_energy){
      khigh = k-1;
      break;
    }
  }

  if(khigh==MAX_ENERGY_BIN) cgmTerminateCode("continuum energy bins too large",khigh);

  n->ntotal = khigh+2;      // number of bins, incl. discrete region

  /*** determine the mid-point energies in each bin.
       the sequence starts from the highest toward zero,
       so that n.excitation[0] = n.max_energy */

  for(int k=0 ; k<n->ntotal ; k++){
#ifdef  HAVE_PRIVATE_ENERGY_GRID
    n->excitation[k] = n->max_energy - custom_energy_grid[k];
#else
    n->excitation[k] = n->max_energy - binenergy(n->binwidth,k) + ENERGY_BIN*0.5;
#endif
  }
  n->excitation[n->ntotal-1] = 0.0;

  /*** find a bin just above discrete level */
  n->ncont  = 0;            // number of bins in continuum only
  if(elmax == 0.0) n->ncont = n->ntotal-1;
  else{
    if(elmax < n->max_energy){
      for(int k=1 ; k<n->ntotal ; k++){
        if(n->excitation[k] < elmax){
          n->ncont = k -1;
          break;
        }
      }
    }
  }
/*
  for(int k=0 ; k<n->ntotal ; k++){
    cout << k << " " << n->excitation[k] << " " << n->binwidth[k]
         << " "<< n->max_energy - n->excitation[k] << endl;
  }
  cout << "N: " << n->ntotal <<" " << n->ncont << " " << n->ndisc
       << " "<< elmax << " "<< khigh << endl;
*/
  return(n->ncont);
}


/***********************************************************/
/*     Prepare Level Density Parameters                    */
/***********************************************************/
void statSetupLevelDensityParameter(Nucleus *n, LevelDensity *ldp)
{
  /*** read in level density parameters */
//  kckDataRead(&n->za, ldp);
  getkcksystdat(&n->za, ldp);
  ldp->a = kckAsymptoticLevelDensity(n->za.getA());


  /*** determine sigma2 from discrete levels */
  ldp->sigma0 = ldLevelSpinCutoff(n->ndisc, n->lev);
  if(ldp->sigma0 == 0.0) ldp->sigma0 = sqrt(kckSpinCutoff(n->za.getA()));

  /*** re-connect ConstTemp and Fermi Gas, because Nlevel can be changed */
  ldTGconnect((double)n->za.getA(), n->lev[n->ndisc-1].energy, (double) n->ndisc, ldp);

  /*** if no connection, use systematics */
  if(ldp->match_energy == 0.0){
    ldp->temperature = kckTemperature(n->za.getA(),ldp->shell_correct);
		
    /*** fluctuate temperature */
    //    ldp->temperature = normaldist(ldp->temperature,0.15);

    ldp->E0 = kckE0(n->za.getA(),ldp->pairing_energy,ldp->shell_correct);
    ldTextrapolate((double)n->za.getA(), n->lev[n->ndisc-1].energy, ldp);
  }

  statSetupLevelDensity(n,ldp);
}


/***********************************************************/
/*     Store Level Densities in the Continuum Array        */
/***********************************************************/
void statSetupLevelDensity(Nucleus *n, LevelDensity *ldp)
{
  for(int k=0 ; k<=n->ncont ; k++){
    double ex = n->excitation[k];
    double r  = ldLevelDensity(ex,(double)n->za.getA(),ldp);

    for(int j=0 ; j<MAX_J ; j++){
      double sd = ldSpinDistribution(j+halfint(n->lev[0].spin),ldp->spin_cutoff,ldp->sigma0);
      n->density[k][j].even = r*sd*ldParityDistribution( 1);
      n->density[k][j].odd  = r*sd*ldParityDistribution(-1);
    }
  }

/*  cout << n->ncont << endl;
  for(int k=0 ; k <= n->ncont ; k++){
  cout << n->excitation[k] << " ";
    for(int j=0 ; j<8 ; j++)  cout << " " << n->density[k][j].even+n->density[k][j].odd;
    cout << endl;
  }
 */

}


/***********************************************************/
/*     Generate Fictitious Levels                          */
/***********************************************************/
void statGenerateLevels(int nlev, Nucleus *n, LevelDensity *ldp)
{
  double *elev;
  elev = new double [MAX_LEVELS];

  /*** determine the highest energy for the discrete level region
       by integrating level densities. 
       max set to 1 MeV, min set to 5 keV */
  double de = 0.001, elmax = 1.0, elmin = 0.005;
  int    km = (int)(elmax/de);

  double s = 0.0;
  for(int k=1 ; k<=km ; k++){
    double ex = k*de;
    if(ex > elmax) break;
    s += ldLevelDensity(ex,(double)n->za.getA(),ldp)*de;
    if(s > (double)nlev+1){ // including ground state
      elmax = ex;
      break;
    }
  }

  /*** distribute N-levels randomly in [elmin,elmax] */
  n->ndisc = nlev + 1;

  for(int i=0 ; i<nlev ; i++)
    elev[i] = (i==nlev) ? elmax : rng_cgm() * (elmax-elmin) + elmin;
  sort(elev,elev+nlev);

  for(int i=1 ; i<=nlev ; i++){
    n->lev[i].energy   = elev[i-1];
    n->lev[i].spin     = -1.0;
    n->lev[i].parity   =  0;
    n->lev[i].halflife =  0.0;
    n->lev[i].ngamma   =  1;
    
    for(int j=0 ; j<n->lev[i].ngamma ; j++){
      if(j==0){
        n->lev[i].fstate[j] = 0;
        n->lev[i].branch[j] = 1.0;
        n->lev[i].gratio[j] = 1.0;
      }
      else{
        n->lev[i].fstate.push_back(0);
        n->lev[i].branch.push_back(1.0);
        n->lev[i].gratio.push_back(1.0);
      }
    }
  }
  statFixDiscreteLevels(n);

  /*** fix number of continuum bins */
  n->ncont = 0;
  if(elmax < n->max_energy){
    for(int k=1 ; k<n->ntotal ; k++){
      if(n->excitation[k] < elmax){
        n->ncont = k -1;
        break;
      }
    }
  }
//  for(int i=0 ; i<n->ndisc ; i++)
//    cout << n->lev[i].energy <<" "<< n->lev[i].spin << " " << n->lev[i].parity << endl;

  delete [] elev;
}


/**********************************************************/
/*      Fix Discrete Level Data If Not Assigned           */
/**********************************************************/
void statFixDiscreteLevels(Nucleus *n)
{  
  /*** quick scan if data are incomplete */
  bool needfix = false;
  for(int i=0 ; i<n->ndisc ; i++){
     if( (n->lev[i].spin < 0.0) || (n->lev[i].parity == 0) ) needfix = true;
  }
  if(!needfix) return;

  for(int i=0 ; i<n->ndisc ; i++){

    /*** re-assign spin, sampling from spin and parity distributions */
    if(n->lev[i].spin < 0.0){
      double s0 = sqrt(kckSpinCutoff(n->za.getA()));
      double r = rng_cgm();
      double s = 0.0;
      double c = sqrt(PI2)*s0*s0*s0;
      for(int j=0 ; j<MAX_J ; j++){
        double x = j+halfint(n->lev[0].spin);
        s += (x+0.5)/c * exp(-(x+0.5)*(x+0.5)/(2*s0*s0)) * (2*x+1.0);
        if(s>r){
          n->lev[i].spin = x;
          break;
        }
      }
    }
    
    /*** assume even distribution for parity */
    if(n->lev[i].parity == 0){
      double r = rng_cgm();
      n->lev[i].parity = (r>0.5) ? 1 : -1;
    }
  }

  return;
}


/***********************************************************/
/*     Get GDR Parameters                                  */
/***********************************************************/
void statSetupGdrParameter(Nucleus *n, double a, GDR *gdr, double beta)
{
  int  m = 0;

  /*** clear GDR parameters */
  for(int i=0 ; i<MAX_GDR ; i++) gdr[i].clear();

  /*** use E1 systematics */
  if(beta == 0.0){
    gdrDArigo((double)n->za.getA(),beta,&gdr[m++]);
  }else{
    gdrDArigoDoubleHump0((double)n->za.getA(),beta,&gdr[m++]);
    gdrDArigoDoubleHump1((double)n->za.getA(),beta,&gdr[m++]);
  }

  /*** M1 and E2 */
  gdrM1((double)n->za.getA(),&gdr[m++]);
  gdrE2((double)n->za.getZ(),(double)n->za.getA(),&gdr[m++]);

  /*** renormalize M1 photo cross section */
  gdrM1norm((double)n->za.getA(),a,gdr);

  /*** scissors mode (M1) or pygmy resonance (E1), if needed */
//gdr[m++].setGDR("M1",2.0,0.6,1.2);
//gdr[m++].setGDR("E1",2.0,0.6,1.2);
  
  /*** save parameters with in gtrans.cpp scope */
  gdrParameterSave(m,gdr);
}


/***********************************************************/
/*      Count Number of Possible Neutrons to Be Emitted    */
/***********************************************************/
int statTotalNeutrons(double emax, ZAnumber *cnZA)
{
  int n = 0;
  int z = cnZA->getZ();
  int a = cnZA->getA();

  for(int i=0 ; i<10 ; i++){
    double sn = mass_excess(z,a-1) + ENEUTRON - mass_excess(z,a);
//  cout << "# Etotal: " << emax << "  Sn: " << sn << endl;
    emax -= sn;
    if(emax <= 0.0) break;
    n++;
    a--;
  }

  if(n >= MAX_COMPOUND) n = MAX_COMPOUND-1;

//cout << "# number of neutrons "<< n << endl;

  return n;
}


/***********************************************************/
/*      Get an Excitation Energy of the Bin                */
/***********************************************************/
double binenergy(double *d, int x){
  double e = 0.0;
  if(x == 0) e = d[x]*0.5;
  else{
    e = d[0]*0.5;
    for(int k=1 ; k<=x ; k++) e += d[k];
  }
  return (e);
}


/***********************************************************/
/*      Normal Distribtion Random Numbers                  */
/***********************************************************/
double normaldist(double m, double s)
{
  double a,b,r;
  do{
    a = 2.0*rng_cgm()-1.0;
    b = 2.0*rng_cgm()-1.0;
    r  = a*a+b*b;
  }while (r>=1.0 || r==0.0);

  b = a*sqrt(-2.0*log(r)/r);

  return(s*b+m);
}
