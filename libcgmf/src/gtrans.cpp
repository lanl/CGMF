/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file gtrans.cpp

  \brief Gamma-ray transmission coefficients

*/

#include <iostream>
#include <cmath>

using namespace std;

#include "physics.h"
#include "cgm.h"

static double  gdrStandardLorentzian    (int, double, GDR *);
static double  gdrGeneralizedLorentzian (int, double, GDR *, double, double);
static inline double  gdrProfileLorentzian     (double, double, double);
static inline double  gdrTempDependGamma       (double, double, double, double);

//static double gStrengthFunction = 1.0;
static double gStrengthFunction = 3.14035; // determined in the FP region

static GDR gdr[MAX_GDR];
static int ngdr = 0;

/**********************************************************/
/*      Normalization Factor for M1 Transmission          */
/**********************************************************/
void gdrM1norm(double a, double ldpa, GDR *g)
{
  double fe1 = 0.0, fm1 = 0.0;

  /*** scan all E1 GDR */
  for(int i=0 ; i<MAX_GDR ; i++){
    if(g[i].getXL() == "E1"){
      fe1 += gdrGeneralizedLorentzian(E1_MULTIPOL,7.0,&g[i],ldpa,0.0);
    }
  }

  /*** look for M1, should be one */
  int m = 0;
  for(int i=0 ; i<MAX_GDR ; i++){
    if(g[i].getXL() == "M1"){
      fm1 += gdrStandardLorentzian(M1_MULTIPOL,7.0,&g[i]);
      m = i;
      break;
    }
  }
  g[m].setSigma(17.01*pow(a,-0.878)*fe1/fm1);
}


/**********************************************************/
/*      Set GDR Parameters in Scope Here                  */
/**********************************************************/
void gdrParameterSave(int nmp, GDR *gdr0)
{
  for(int i=0 ; i<nmp ; i++){
    if(i<MAX_GDR) gdr[i] = gdr0[i];
  }
  ngdr = nmp;
}


/**********************************************************/
/*      Gamma-ray Transmission                            */
/**********************************************************/
double gdrGammaTransmission(GammaProfile c, int l, double eg, GammaMultipol gm, double a, double ex)
{
  string xl = "  ";

  switch(gm){
  case E1: xl = "E1"; break;
  case M1: xl = "M1"; break;
  case E2: xl = "E2"; break;
  default:            break;
  }

  double tg = 0.0;
  for(int i=0 ; i<ngdr ; i++){
    if(gdr[i].getXL() == xl){
      if(c==SL){
        tg += gdrStandardLorentzian(l,eg,&gdr[i]);
      }
      else if(c==GL){
        tg += gdrGeneralizedLorentzian(l,eg,&gdr[i],a,ex);
      }
    }
  }
  tg *= PI2*pow(eg,2*l+1.0)*gStrengthFunction;

  return (tg);
}


/**********************************************************/
/*      Brink-Axel Standard Lorentzian Function           */
/**********************************************************/
double gdrStandardLorentzian(int l, double eg, GDR *gam)
{
  double c,f;

  c = 26.0e-08/(2*l+1.0);
  f = c * gam->getSigma() * gam->getWidth() * pow(eg,2.0-2.0*l)
        * gdrProfileLorentzian(eg,gam->getEnergy(),gam->getWidth());

   return(f);
}


/**********************************************************/
/*      Kopecky-Uhl Generalized Lorentzian Function       */
/**********************************************************/
double gdrGeneralizedLorentzian(int l, double eg, GDR *gam, double a, double ex)
{
  double c,g,f,t;

  c = 26.0e-08/(2*l+1.0);
  t = (ex<0.0) ? 0.0 : sqrt(ex/a);
  g = (l==1) ? gdrTempDependGamma(eg,gam->getEnergy(),gam->getWidth(),t) 
             : gam->getWidth();
  f = c*gam->getSigma() * gam->getWidth() * pow(eg,2.0-2.0*l)
      *( gdrProfileLorentzian(eg,gam->getEnergy(),g) 
      +( (l==1) ? 0.7*gdrTempDependGamma(0.0,gam->getEnergy(),gam->getWidth(),t)
                  /(gam->getEnergy() * gam->getEnergy() * gam->getEnergy())
                : 0.0 )  );
/*
   printf(" % 11.4e % 11.4e % 11.4e % 11.4e % 11.4e % 11.4e\n",
          eg,gdrProfileLorentzian(eg,gam->getEnergy(),g),t,g,f,PI2*pow(eg,2*l+1.0)*f);
 */
  return(f);
}


/**********************************************************/
/*      Temperature dependent Gamma, by Kadmenskij        */
/**********************************************************/
double gdrTempDependGamma(double eg, double er, double w, double t)
{
  if(er==0.0) return(0.0);
  return( w*(eg*eg + 4.0*PI*PI*t*t)/(er*er) );
}


/**********************************************************/
/*      Lorentzian Profile                                */
/**********************************************************/
double gdrProfileLorentzian(double eg, double er, double w)
{
  if(er==eg && eg*w==0.0) return(0.0);
  return (eg*w /( (er*er-eg*eg)*(er*er-eg*eg) + eg*eg*w*w )); 
}
