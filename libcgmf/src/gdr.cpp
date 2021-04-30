/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file gdr.cpp

  \brief Global GDR parameter library

*/

#include <iostream>
#include <cmath>

using namespace std;

#include "gdr.h"


/**********************************************************/
/*      E1                                                */
/**********************************************************/
void gdrDArigo(double mass, double beta, GDR *g)
{
  double energy = (49.336+7.34*beta)*pow(mass,-0.2409);
  double width  = 0.3 * energy;
  double sigma  = 10.6 * mass / width;
  string XL     = "E1";
  g->setGDR(XL,energy,width,sigma);
}

void gdrDArigoDoubleHump0(double mass, double beta, GDR *g)
{
  double energy = 50.0*pow(mass,-0.232) * exp(-0.946*beta);
  double width  = (0.282-0.263*beta) * energy;
  double sigma  = 3.48  * mass / width;
  string XL     = "E1";
  g->setGDR(XL,energy,width,sigma);
}

void gdrDArigoDoubleHump1(double mass, double beta, GDR *g)
{
  double energy = 50.0*pow(mass,-0.232);
  double width  = (0.35 -0.14 *beta) * energy;
  double sigma  = 8.26  * mass / width;
  string XL     = "E1";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      M1                                                */
/**********************************************************/
void gdrM1(double mass, GDR *g)
{
  double energy = 41.0/pow(mass,1/3.0);
  double width  = 4.0;
  double sigma  = 1.0;
  string XL     = "M1";
  g->setGDR(XL,energy,width,sigma);
}


/**********************************************************/
/*      E2                                                */
/**********************************************************/
void gdrE2(double z, double mass, GDR *g)
{
  double a3 = 1./pow(mass,1/3.0);
  double energy = 63.0*a3;
  double width  = 6.11 - 0.012*mass;
  double sigma  = 1.5e-4 * z * z * energy * energy*a3 / width;
  string XL     = "E2";
  g->setGDR(XL,energy,width,sigma);
}
