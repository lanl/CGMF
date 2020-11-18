/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file excinterface.cpp

  \brief Exciton model wrapper

*/

#include <iostream>

using namespace std;

#include "physics.h"
#include "cgm.h"
#include "masstable.h"
#include "excinterface.h"
#include "terminate.h"


void excitonInterface(const int targZ, const int targA, const double elab, double **spc)
{
  Transmission  *tc[MAX_CHANNEL];
  Pdata         pdt[MAX_CHANNEL];
  System        sys;

  try{
    for(int j=1 ; j<MAX_CHANNEL ; j++){
      tc[j] = new Transmission [MAX_ENERGY_BIN];
      for(int k=0 ; k<MAX_ENERGY_BIN ; k++) tc[j][k].tran = new double [3*MAX_J];
    }
  }
  catch(bad_alloc){
    cgmTerminateCode("memory allocation error in excinterface");
  }

  for(int k=0 ; k<MAX_ENERGY_BIN ; k++) spc[gammaray][k] = spc[neutron][k] = 0.0;


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

  /*** copy input variables */
  sys.lab_energy = elab;
  sys.target.setZA(targZ,targA);
  sys.compound.setZA(targZ,targA+1);
  sys.incident = pdt[neutron];

  ncl[0].za = sys.compound;
  ncl[1].za = sys.target;

  /*** calculate maximum energy of the system */
  double m0 = mass_excess(sys.compound.getZ(),sys.compound.getA());
  double m1 = mass_excess(sys.target.getZ()  ,sys.target.getA()  );
  double m2 = sys.target.getA() + m1/AMUNIT;
  double sn = m1 + ENEUTRON - m0;

  sys.cms_energy = sys.lab_energy * m2 / (m2 + NEUTRONMASS);
  sys.ex_total   = sn + sys.cms_energy;
  ncl[0].max_energy = sys.ex_total;

  /*** initialize nucleus and particle data */
  int nemit = 1;
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

  for(int i=1 ; i<=nemit ; i++){
    double sepa = ncl[i].mass_excess + pdt[neutron].mass_excess - ncl[i-1].mass_excess;
    ncl[i-1].cdt[neutron].binding_energy = sepa;

    if(sepa < ncl[i-1].max_energy){
      ncl[i-1].cdt[neutron].status = true;
      ncl[i-1].cdt[neutron].next   = i;
      ncl[i  ].max_energy          = ncl[i-1].max_energy - sepa;
    }
  }

  ncl[0].cdt[0].next   = 0;
  ncl[0].cdt[0].status = true;


  /*** set up energy bins */ 
  for(int i=0 ; i<=nemit ; i++){
    statSetupEnergyBin(&ncl[i]);
  }

  /*** calculate inverse reaction cross sections */
  statStoreContinuumTransmission(0,0,&pdt[neutron],tc[neutron]);

  /*** exciton model */
  preqExcitonModel(ncl, &sys, pdt, tc, spc);

  /*** normalize pre-equilibrium spectrum */
  double sum = 0.0;
  for(int k=0 ; k<ncl[1].ntotal ; k++) sum += spc[neutron][k];
  if(sum > 0.0){
      for(int k=0 ; k<ncl[1].ntotal ; k++) {spc[neutron][k] /= sum;}
  }

  for(int j=1 ; j<MAX_CHANNEL ; j++){
    for(int k=0 ; k<MAX_ENERGY_BIN ; k++) delete [] tc[j][k].tran;
    delete [] tc[j];
  }
}



