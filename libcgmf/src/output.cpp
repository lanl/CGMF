/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file output.cpp

  \brief Main data output

*/

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cmath>

using namespace std;

#include "physics.h"
#include "cgm.h"
#include "global.h"
//#include "global_var.h"

static const double eps = 1.0e-99;

static inline double lowfilter (double x){if(fabs(x)<eps) return 0.0; else return(x);}

static inline double cgmGaussianBroadening (int, int, double, double *);

/***********************************************************/
/*      Set Gamma Spectra to pass back to mcnpx            */
/***********************************************************/
void cgmGetSpectra(double de, double **spc)
{
  int    cm = 2, k0 = 0;

  mcl_nlines=0;

  for(int k=0 ; k<=MAX_ENERGY_BIN-1 ; k++){
    double emin = (k>0)  ? ((double)k-0.5)*de : 0;
    double emax = ((double)k+0.5)*de;
    double e    = (emin+emax)/2.0;
    double d    =  emax-emin;

      for(int i=0 ; i<cm ; i++){
	if (lowfilter(spc[i][k]) > 0.0) {
          if(mcl_nlines < 30 ) mcl_glines[mcl_nlines] = (emin+emax)/2.0;
          mcl_nlines++;
        }
      }
  }
}

/***********************************************************/
/*      Print Gamma, Neutron, Electron, Neutrino Spectra   */
/***********************************************************/
void cgmPrintSpectra(bool beta, double de, double **spc)
{
  double s[4], u[4], x[4];
  int    cm = 2, k0 = 0;

  if(beta) cm = 4; // include electron and neutrino channels

  k0 = cgmZeroCut(cm,spc);

  /*** integration for nomalization */
  for(int i=0 ; i<cm ; i++){
    s[i] = 0.0;
    for(int k=0 ; k<=k0 ; k++){
      double emin = (k>0)  ? ((double)k-0.5)*de : 0;
      double emax = ((double)k+0.5)*de;
      double d    =  emax-emin;
      s[i]       += spc[i][k];
      spc[i][k]  /= d;
    }
    u[i] = (s[i]>0.0) ? 1.0/s[i] : 0.0;
  }

  cout << "#      Emin        Emax        Ecal        "
       << "  Spectra/Decay/MeV          " << endl;
  cout << "#                                        " << "  Gamma-Ray   Neutron   ";
  if(beta) cout << "  Electron    Neutrino  ";
  cout << endl;

  for(int i=0 ; i<cm ; i++) x[i] = 0.0;
  for(int k=0 ; k<=k0 ; k++){
    double emin = (k>0)  ? ((double)k-0.5)*de : 0;
    double emax = ((double)k+0.5)*de;
    double e    = (emin+emax)/2.0;
    double d    =  emax-emin;

    for(int i=0 ; i<cm ; i++) x[i] += spc[i][k]*e*d;

    cout << setw(5) << k;
    cout << setprecision(4) << setiosflags(ios::scientific);
    cout << setw(12) << emin
         << setw(12) << emax
         << setw(12) << e;

    if(ctl.print_broadened){
      cout << setw(12) << lowfilter( cgmGaussianBroadening(k0,k,de,spc[0]) );
      for(int i=1 ; i<cm ; i++) cout << setw(12) << lowfilter(spc[i][k]);
    }else{
      for(int i=0 ; i<cm ; i++) cout << setw(12) << lowfilter(spc[i][k]);
    }
    cout << endl;
  }

  cout << "#----------------------------------------------------------------" << endl;
  cout << setprecision(4) << setiosflags(ios::scientific);

  cout << "# E (total)                              ";
  for(int i=0 ; i<cm ; i++) cout << setw(12) << lowfilter(x[i]);
  cout << endl;

  cout << "# E (average)                            ";
  for(int i=0 ; i<cm ; i++) cout << setw(12) << lowfilter(x[i]*u[i]);
  cout << endl;

  cout << "# Multiplicity                           ";
  for(int i=0 ; i<cm ; i++) cout << setw(12) << lowfilter(s[i]);
  cout << endl;

  cout << endl;
  cout << endl;
  /*
  double etot = 0.0;
  for(int i=0 ; i<cm ; i++) etot += x[i];
  cout << "# Total Energy Release       " << setw(12) << etot << endl;
  */
}


/***********************************************************/
/*      Gaussian Broadening                                */
/***********************************************************/
inline double cgmGaussianBroadening(int k0, int k, double de, double *spec)
{
  const double d = 0.05; // energy resolution 50keV
  const double c = 1.0/(sqrt(PI2)*d);

  double r = 0.0;
  for(int j=0 ; j<=k0 ; j++){
    double e  = (j-k)*de;
    r += c * spec[j] * exp(-e*e/(2*d*d))*de;
  }
  return(r);
}


/***********************************************************/
/*      Eliminate Zeros at High Energies                   */
/***********************************************************/
int cgmZeroCut(int cm, double **spc)
{
  int k0 = 0;

  for(int k=MAX_ENERGY_BIN-1 ; k>0 ; k--){
    bool nonzero = false;
    for(int c=0 ; c<cm ; c++) if(spc[c][k] > 0.0) nonzero = true;
    if(nonzero){
      k0 = k+1;
      break;
    }
  }

  return(k0);
}


/***********************************************************/
/*      Gamma-ray Spectrum for Custom Energy Grid          */
/***********************************************************/
void cgmPrintCustomGrid(int ns, double *es, double *de, double **spc)
{
  double s[2], emin, emax;

  cout << setprecision(4) << setiosflags(ios::scientific);

  cout << "#      Emin        Emax        Ecal        "
       << "Gamma-Ray   "
       << "Neutron     " << endl;

  s[0] = s[1] = 0.0;
  emin = 0.0;
  emax = de[0]*0.5;
  for(int k=0 ; k<ns ; k++){
    s[0] += spc[0][k];
    s[1] += spc[1][k];

    spc[0][k] /= de[k];
    spc[1][k] /= de[k];

    cout << setw(5) << k;
    cout << setprecision(4) << setiosflags(ios::scientific);
    cout << setw(12) << emin
         << setw(12) << emax
         << setw(12) << es[k];
    cout << setw(12) << lowfilter(spc[0][k])
         << setw(12) << lowfilter(spc[1][k]) << endl;

    emin = emax;
    emax += de[k+1];
  }
  cout << "#----------------------------------------------------------------" << endl;
  cout << setprecision(4) << setiosflags(ios::scientific);

  cout << "# Sum                                    ";
  for(int i=0 ; i<2 ; i++) cout << setw(12) << lowfilter(s[i]);
//cout << setw(12) << lowfilter(s[0]+s[1]) << endl;
  cout << endl;

  cout << endl;
  cout << endl;
}


/***********************************************************/
/*      Output Population in the Residual Nucleus          */
/***********************************************************/
void cgmPrintPopulation(Nucleus *n0, Nucleus *n1)
{
  cout << "# Gam CONT  " << setw(12) << n0->ncont << setw(12) << n0->jmax+1 << endl;
  cout << "# ";
  for(int j=0 ; j<=n0->jmax ; j++){
    cout << setw(5) << j << "  (+)  ";
    cout << setw(5) << j << "  (-)  ";
  }
  cout << endl;

  for(int k=0 ; k<n0->ncont ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(8) << setw(12) << n0->excitation[k];
    cout << setiosflags(ios::scientific) << setprecision(5);

    if(k == 0){
      /*** because initial population is stored at k=0, we skip it */
      for(int j=0 ; j<=n0->jmax ; j++) cout << setw(12) << 0.0 << setw(12) << 0.0;
    }
    else{
      for(int j=0 ; j<=n0->jmax ; j++) {
        cout << setw(12) << n0->pop[k][j].even << setw(12) << n0->pop[k][j].odd;
      }
    }
    cout << endl;
  }

  cout << "# Gam DISC  " << setw(12) << n0->ndisc << endl;
  for(int k=0 ; k<n0->ndisc ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(8) << setw(12) << n0->lev[k].energy;
    cout << setiosflags(ios::scientific) << setprecision(5);
    cout << setw(12) << n0->lpop[k] << endl;
  }
  cout << endl;

  cout << "# Neu CONT  " << setw(12) << n1->ncont << setw(12) << n1->jmax+1 << endl;
  cout << "# ";
  for(int j=0 ; j<=n0->jmax ; j++){
    cout << setw(5) << j << "  (+)  ";
    cout << setw(5) << j << "  (-)  ";
  }
  cout << endl;

  for(int k=0 ; k<n1->ncont ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(8) << setw(12) << n1->excitation[k];
    cout << setiosflags(ios::scientific) << setprecision(5);

    for(int j=0 ; j<=n1->jmax ; j++) {
      cout << setw(12) << n1->pop[k][j].even << setw(12) << n1->pop[k][j].odd;
    }
    cout << endl;
  }

  cout << "# Neu DISC  " << setw(12) << n1->ndisc << endl;
  for(int k=0 ; k<n1->ndisc ; k++){
    cout.setf(ios::fixed, ios::floatfield);
    cout << setprecision(8) << setw(12) << n1->lev[k].energy;
    cout << setiosflags(ios::scientific) << setprecision(5);
    cout << setw(12) << n1->lpop[k] << endl;
  }
  cout << endl;
}
