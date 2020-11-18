//
//  cgmf.cpp
//

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <cstdlib>

#include "cgmfEvents.h"
#include "rngcgm.h"
#include "config.h"
#include "config-ff.h"

#include "cgmf_config.h"

#ifdef MPIRUN
#include <mpi.h>
#endif

using namespace std;

void readUserInput      (int, char *[], int);
void recordEvent        (cgmfEvent *);
void printEventToFile   (FILE *fp, cgmfEvent *, double);
void printSummaryEvents (cgmfEvent *);

// Default user input values
string path = "";
double incidentEnergy=-1;
int ZAIDt=0;
int nevents=0;
int startingEvent=1;
string outfilename="";
double timeCoincidenceWindow=1e-8;

int    sumALF=0, sumAHF=0, sumZLF=0, sumZHF=0;
int    sumNuLF=0, sumNuHF=0, sumNuPre=0, sumNuTot=0;
int    sumNugLF=0, sumNugHF=0, sumNugTot=0;
double sumKELF=0.0, sumKEHF=0.0, sumTKE=0.0;
double sumXELF=0.0, sumXEHF=0.0, sumTXE=0.0;
double sumJLF=0.0, sumJHF=0.0, sumJ=0.0;

double sumEcmLF=0.0, sumEcmHF=0.0;
double sumElabLF=0.0, sumElabHF=0.0;
double sumEnPre=0.0;
double sumEcm=0.0, sumElab=0.0;
double sumEgLF=0.0, sumEgHF=0.0, sumEg=0.0;

//////////////////////////////////////////////////////////////////////////////////////
//                                   MAIN DRIVER 
//////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[]) {

  int ip=0;

  cgmfEvent*  event = 0;
  cgmfYields* yields;

#ifdef MPIRUN
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
#endif

  readUserInput(argc, argv, ip); // check the validity of the user input

  // initialize random number generator
  UniformRNG rng(1);
	
  if (nevents<0) {
    if(ip==0){
#ifdef MPIRUN
      cout << "Not implemented: Fission fragment yields calculation (negative -nevents option) will run with a single MPI rank.\n";
#endif
      rng.set_seed(startingEvent);
      set_rng(rng);
      yields = new cgmfYields (ZAIDt, incidentEnergy, -nevents, outfilename);
      printf("\n/// CGMF-generated scission fragment yields Y(Z,A,KE,U,J,Pi,px,py,pz) saved in file %s ///\n", outfilename.c_str());
    }
  } else {
    FILE *fp = fopen(&outfilename[0],"w");
    if(ip==0) fprintf(fp, "# %5i %g %g\n", ZAIDt, incidentEnergy,timeCoincidenceWindow);
    for (int i=0; i<nevents; i++) {
      rng.set_seed(i+ip*nevents+startingEvent);
      set_rng(rng);
//      if (i>9 and i%(nevents/10)==0 and ip==0) cout << float(i)/float(nevents)*100.0 << "%\n";
      if (event != 0) delete event;
      event = new cgmfEvent(ZAIDt, incidentEnergy, 0.0, timeCoincidenceWindow);
      recordEvent (event);
      printEventToFile (fp, event, timeCoincidenceWindow);
    }
    fclose(fp);
    printSummaryEvents (event);
  }

  cgmf_cleanup();
  if (event != 0) delete event;
#ifdef MPIRUN
  MPI_Finalize();
#endif

  return 0;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////
//                                        THE END
//////////////////////////////////////////////////////////////////////////////////////////////



/*!
------------------------------------------------------------------------------
\brief Record a CGMF event into arrays to provide summary of N events
------------------------------------------------------------------------------
*/
void recordEvent (cgmfEvent *event) {

  static int i=0;

  int nuLF, nuHF, nuPre, nuTot;
  int nugLF, nugHF, nugTot;

  // Light fragment (Z,A)
  sumALF += event->getLightFragmentMass();
  sumZLF += event->getLightFragmentCharge();

  // Heavy fragment (Z,A)
  sumAHF += event->getHeavyFragmentMass();
  sumZHF += event->getHeavyFragmentCharge();

  // Kinetic energies (MeV)
  sumKELF += event->getLightFragmentKineticEnergy();
  sumKEHF += event->getHeavyFragmentKineticEnergy();
  sumTKE   = sumKELF+sumKEHF;

  // Excitation energies (MeV)
  sumXELF += event->getLightFragmentExcitationEnergy();
  sumXEHF += event->getHeavyFragmentExcitationEnergy();
  sumTXE   = sumXELF+sumXEHF;

  // Spin (hbar)
  sumJLF += event->getLightFragmentSpin();
  sumJHF += event->getHeavyFragmentSpin();
  sumJ    = (sumJLF+sumJHF)/2;

  // Neutron multiplicities (n/f)
  nuLF  = event->getLightFragmentNeutronNu();
  nuHF  = event->getHeavyFragmentNeutronNu();
  nuPre = event->getPreFissionNeutronNu();

  sumNuLF  += nuLF;
  sumNuHF  += nuHF;
  sumNuPre += nuPre;
  sumNuTot  = sumNuLF + sumNuHF + sumNuPre;

  // Neutron energies emitted from light fragment, in center-of-mass and lab frame (MeV)
  for (int j=0; j<nuLF; j++) { 
    sumEcmLF  += event->getCmNeutronEnergy(j);
    sumElabLF += event->getNeutronEnergy(j);
  }

  // Neutron energies emitted from heavy fragment, in center-of-mass and lab frame (MeV)
  for (int j=nuLF; j<nuLF+nuHF; j++) {
    sumEcmHF  += event->getCmNeutronEnergy(j);
    sumElabHF += event->getNeutronEnergy(j);
  }

  // Pre-fission neutron kinetic energies (MeV)
  for (int j=0; j<nuPre; j++) sumEnPre += event->getPreFissionNeutronEnergy(j);

  sumEcm  = sumEcmLF + sumEcmHF + sumEnPre;
  sumElab = sumElabLF + sumElabHF + sumEnPre;

  // Photon multiplicities (g/f)
  nugLF  = event->getLightFragmentPhotonNu();
  nugHF  = event->getHeavyFragmentPhotonNu();
  nugTot = nugLF+nugHF;

  sumNugLF += nugLF;
  sumNugHF += nugHF;
  sumNugTot = sumNugLF + sumNugHF;

  // Gamma energies (MeV)
  for (int j=0; j<nugLF; j++)           sumEgLF += event->getPhotonEnergy(j); // Light fragment
  for (int j=nugLF; j<nugLF+nugHF; j++) sumEgHF += event->getPhotonEnergy(j); // Heavy fragment

  sumEg = sumEgLF + sumEgHF;

  i++;

  return;

}

/*!
------------------------------------------------------------------------------
\brief Print a summary of all CGMF events
------------------------------------------------------------------------------
*/
void printSummaryEvents (cgmfEvent* event) {

  int ip=0,np=1;

#ifdef MPIRUN

  MPI_Comm_rank (MPI_COMM_WORLD, &ip);
  MPI_Comm_size (MPI_COMM_WORLD, &np);
  int itmp;
  MPI_Reduce(&sumZLF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumZLF=itmp;
  MPI_Reduce(&sumALF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumALF=itmp;
  MPI_Reduce(&sumZHF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumZHF=itmp;
  MPI_Reduce(&sumAHF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumAHF=itmp;
  double rtmp;
  MPI_Reduce(&sumKELF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumKELF=rtmp;
  MPI_Reduce(&sumKEHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumKEHF=rtmp;
  MPI_Reduce(&sumTKE, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumTKE=rtmp;  
  MPI_Reduce(&sumXELF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumXELF=rtmp;
  MPI_Reduce(&sumXEHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumXEHF=rtmp;
  MPI_Reduce(&sumTXE, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumTXE=rtmp;  
  MPI_Reduce(&sumJLF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumJLF=rtmp;
  MPI_Reduce(&sumJHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumJHF=rtmp;
  MPI_Reduce(&sumJ, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumJ=rtmp;  

  MPI_Reduce(&sumNuLF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNuLF=itmp;
  MPI_Reduce(&sumNuHF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNuHF=itmp;
  MPI_Reduce(&sumNuPre, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNuPre=itmp;
  MPI_Reduce(&sumNuTot, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNuTot=itmp;  
  MPI_Reduce(&sumEcmLF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEcmLF=rtmp;
  MPI_Reduce(&sumEcmHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEcmHF=rtmp;
  MPI_Reduce(&sumEcm, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEcm=rtmp;  
  MPI_Reduce(&sumElabLF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumElabLF=rtmp;
  MPI_Reduce(&sumElabHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumElabHF=rtmp;
  MPI_Reduce(&sumElab, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumElab=rtmp;
  MPI_Reduce(&sumEnPre, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEnPre=rtmp;

  MPI_Reduce(&sumNugLF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNugLF=itmp;
  MPI_Reduce(&sumNugHF, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNugHF=itmp;
  MPI_Reduce(&sumNugTot, &itmp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumNugTot=itmp;  
  MPI_Reduce(&sumEgLF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEgLF=rtmp;
  MPI_Reduce(&sumEgHF, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEgHF=rtmp;
  MPI_Reduce(&sumEg, &rtmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD ) ;
  if(ip==0)
    sumEg=rtmp;

  nevents=np*nevents;
#endif

  if(ip==0){
 
    printf("\n\n //// CGMF Results ////\n\n");

    int ZAIDt=event->getTargetNucleus();
    int Zt=int(ZAIDt/1000.0);
    int At=ZAIDt-1000.0*Zt;

    if (incidentEnergy==0.0) {
      printf("Reaction: spontaneous fission of (%2i,%3i)\n\n", Zt, At);
    } else {
      printf("Reaction: (n,f) on (%2i,%3i) at En=%7.2e MeV\n\n", Zt, At, event->getIncidentEnergy());
    }

    printf("Average Light Fragment (Z,A) = (%.2f,%.2f)\n", double(sumZLF)/nevents, double(sumALF)/nevents);
    printf("Average Heavy Fragment (Z,A) = (%.2f,%.2f)\n\n", double(sumZHF)/nevents, double(sumAHF)/nevents);

    printf("Average Kinetic Energies: LF = %.2f MeV ; HF = %.2f MeV ; <TKE> = %.2f MeV\n",
	   sumKELF/nevents,sumKEHF/nevents,sumTKE/nevents);

    printf("Average Excitation Energies: LF = %.2f MeV ; HF = %.2f MeV ; <TXE> = %.2f MeV\n\n",
	   sumXELF/nevents,sumXEHF/nevents,sumTXE/nevents);

    printf("Average Fragment Spins: <J>_LF = %.2f hbar ; <J>_HF = %.2f hbar ; <J> = %.2f hbar\n",
	   sumJLF/nevents,sumJHF/nevents,sumJ/nevents);

    printf("\n*** Prompt Fission Neutrons ***\n\n");

    printf("Multiplicities (n/f):  <nu>_LF = %.2f ; <nu>_HF = %.2f ; <nu>_prefission = %.2f ; <nu>_tot = %.2f \n",
	   double(sumNuLF)/nevents,double(sumNuHF)/nevents,double(sumNuPre)/nevents,double(sumNuTot)/nevents);

    printf("c-o-m Energies:  <Ecm>_LF = %.2f MeV ; <Ecm>_HF = %.2f MeV ; <Ecm>_prefission = %.2f MeV ; <Ecm>_tot = %.2f MeV\n",
	   sumEcmLF/max(sumNuLF,1),sumEcmHF/max(sumNuHF,1),sumEnPre/max(sumNuPre,1),sumEcm/max(sumNuTot,1));

    printf("Lab. Energies:   <Elab>_LF = %.2f MeV ; <Elab>_HF = %.2f MeV ; <Elab>_prefission = %.2f MeV ; <Elab>_tot = %.2f MeV\n",
	   sumElabLF/max(sumNuLF,1),sumElabHF/max(sumNuHF,1),sumEnPre/max(sumNuPre,1),sumElab/max(sumNuTot,1));

    printf("\n*** Prompt Fission Gammas ***\n\n");

    printf("Multiplicities (g/f):  <nu_g>_LF = %.2f ; <nu_g>_HF = %.2f ; <nu_g>_tot = %.2f \n",
	   double(sumNugLF)/nevents,double(sumNugHF)/nevents,double(sumNugTot)/nevents);

    printf("Gamma Energies:   <Eg>_LF = %.2f MeV ; <Eg>_HF = %.2f MeV ; <Eg>_tot = %.2f MeV\n",
	   sumEgLF/max(sumNugLF,1),sumEgHF/max(sumNugHF,1),sumEg/max(sumNugTot,1));

    printf("\n\n //// THE END ////\n\n");

  }

  return;

}


/*!
------------------------------------------------------------------------------
\brief read user input from the command line, set data directory path, and 
initialize random number generator.
------------------------------------------------------------------------------
*/
void readUserInput (int argc, char *argv[], int ip) {

  int p;

  // read user input from command line
  while ((p=getopt(argc,argv,"e:n:i:f:t:d:s:"))!=-1) {
    switch(p){
      case 'e':
        incidentEnergy=atof(optarg);
        break;
      case 'n':
        nevents=atoi(optarg);
        break;
      case 'i':
        ZAIDt=atoi(optarg);
        break;
      case 'f':
        outfilename=optarg;
        break;
      case 't':
        timeCoincidenceWindow=atof(optarg);
        break;
      case 'd':
        path=optarg;
        break;
      case 's':
        startingEvent=atoi(optarg);
        break;
      default:
        break;
    }
  }

  stringstream sip;
  sip << ip;
  if (outfilename=="") {
  	if (nevents>0) { outfilename="histories.cgmf"; } else { outfilename="yields.cgmf"; }
  }
  outfilename += "."+sip.str();

  if(timeCoincidenceWindow<0.){
    if(timeCoincidenceWindow!=-1.){
      cerr << "The time coincidence window can only be -1 (infinite time coincidence window) or positive" << endl;
      cerr << "Execution terminated" << endl;
#ifdef MPIRUN
      MPI_Finalize();
#endif
      exit(-1);
    }
  }


  // set data directory path
  if (path == "") {
    if (getenv("CGMFDATA") != NULL) {
      path = string(getenv("CGMFDATA"));
    } else if (checkdatapath(INSTALL_DATADIR)) {
      path = INSTALL_DATADIR;
    } else if (checkdatapath(BUILD_DATADIR)) {
      path = BUILD_DATADIR;
    } else {
      cerr << "Cannot find valid CGMFDATA path ... returning" << endl;
      exit(-1);
    }
  }
  if( path.back() != '/' ) path += "/";
  setdatapath(path);


  return;

}

/*!
------------------------------------------------------------------------------
\brief Print a fission event in an output file.
------------------------------------------------------------------------------
*/
void printEventToFile (FILE *fp, cgmfEvent *event, double timeCoincidenceWindow) {

  double eng; // energy
  double diru, dirv, dirw; // directional cosines
  double age; // time of emission of a photon

  // LIGHT FRAGMENT --------------------------

  int nul  = event->getLightFragmentNeutronNu();
  int nugl = event->getLightFragmentPhotonNu();

  fprintf(fp, " %i %i %.3f %.1f %i %.3f %.3f %i %i 0\n", event->getLightFragmentMass(), 
    event->getLightFragmentCharge(), event->getLightFragmentExcitationEnergy(), 
    event->getLightFragmentSpin(), event->getLightFragmentParity(), 
    event->getLightFragmentKineticEnergy(), event->getLightFragmentKineticEnergyPost(), 
    event->getLightFragmentNeutronNu(), event->getLightFragmentPhotonNu());

  fprintf(fp, " %.3f %.3f %.3f %.3f %.3f %.3f\n", event->getLightFragmentPreMomentumX(), 
    event->getLightFragmentPreMomentumY(), event->getLightFragmentPreMomentumZ(), 
    event->getLightFragmentPostMomentumX(), event->getLightFragmentPostMomentumY(), 
    event->getLightFragmentPostMomentumZ());

  // Center-of-Mass neutron (energy, directional cosines) from light fragment
  for (int n1=0; n1<nul; n1++) {
    eng  = event->getCmNeutronEnergy(n1);
    diru = event->getCmNeutronDircosu(n1);
    dirv = event->getCmNeutronDircosv(n1);
    dirw = event->getCmNeutronDircosw(n1);
    fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
  }
  if (nul>0) fprintf(fp, "\n");

  // LAB neutron (energy, directional cosines) from light fragment
  for (int n1=0; n1<nul; n1++) {
    eng  = event->getNeutronEnergy(n1);
    diru = event->getNeutronDircosu(n1);
    dirv = event->getNeutronDircosv(n1);
    dirw = event->getNeutronDircosw(n1);
    fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
  }
  if (nul>0) fprintf(fp, "\n");

  // Photons from light fragment
  for(int n1=0; n1<nugl; n1++) {
    eng  = event->getPhotonEnergy(n1);
    diru = event->getPhotonDircosu(n1);
    dirv = event->getPhotonDircosv(n1);
    dirw = event->getPhotonDircosw(n1);
    if (timeCoincidenceWindow<0) {
      age  = event->getPhotonAge(n1);
      fprintf(fp, "%.3f %.3f %.3f %.3f %g ", diru, dirv, dirw, eng, age);
    } else {
      fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
    }
  }
  if (nugl>0) fprintf(fp, "\n");

  // HEAVY FRAGMENT --------------------------
  int nuh  = event->getHeavyFragmentNeutronNu();
  int nugh = event->getHeavyFragmentPhotonNu();
  int nup  = event->getPreFissionNeutronNu();

  fprintf(fp, " %i %i %.3f %.1f %i %.3f %.3f %i %i %i\n", event->getHeavyFragmentMass(),
    event->getHeavyFragmentCharge(), event->getHeavyFragmentExcitationEnergy(), 
    event->getHeavyFragmentSpin(), event->getHeavyFragmentParity(), 
    event->getHeavyFragmentKineticEnergy(), event->getHeavyFragmentKineticEnergyPost(), 
    event->getHeavyFragmentNeutronNu(), event->getHeavyFragmentPhotonNu(), nup);

  fprintf(fp, " %.3f %.3f %.3f %.3f %.3f %.3f\n", event->getHeavyFragmentPreMomentumX(), 
    event->getHeavyFragmentPreMomentumY(), event->getHeavyFragmentPreMomentumZ(), 
    event->getHeavyFragmentPostMomentumX(), event->getHeavyFragmentPostMomentumY(), 
    event->getHeavyFragmentPostMomentumZ());

  // Center-of-Mass neutron (energy, directional cosines) from heavy fragment
  for (int n1=nul; n1<nul+nuh; n1++) {
    eng  = event->getCmNeutronEnergy(n1);
    diru = event->getCmNeutronDircosu(n1);
    dirv = event->getCmNeutronDircosv(n1);
    dirw = event->getCmNeutronDircosw(n1);
    fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
  }
  if (nuh>0) fprintf(fp, "\n");

  // LAB neutron (energy, directional cosines) from heavy fragment
  for (int n1=nul; n1<nul+nuh; n1++) {
    eng  = event->getNeutronEnergy(n1);
    diru = event->getNeutronDircosu(n1);
    dirv = event->getNeutronDircosv(n1);
    dirw = event->getNeutronDircosw(n1);
    fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
  }
  if (nuh>0) fprintf(fp, "\n");

  // Photons from heavy fragment
  for(int n1=nugl; n1<nugl+nugh; n1++) {
    eng  = event->getPhotonEnergy(n1);
    diru = event->getPhotonDircosu(n1);
    dirv = event->getPhotonDircosv(n1);
    dirw = event->getPhotonDircosw(n1);
    if (timeCoincidenceWindow<0) {
      age  = event->getPhotonAge(n1);
      fprintf(fp, "%.3f %.3f %.3f %.3f %g ", diru, dirv, dirw, eng, age);
    } else {
      fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
    }
  }
  if (nugh>0) fprintf(fp, "\n");

  if (event->getPreFissionNeutronNu()>0) {
    for (int i=0; i<event->getPreFissionNeutronNu(); i++) {
      eng  = event->getPreFissionNeutronEnergy(i);
      diru = event->getPreFissionNeutronDircosu(i);
      dirv = event->getPreFissionNeutronDircosv(i);
      dirw = event->getPreFissionNeutronDircosw(i);
      fprintf(fp, "%.3f %.3f %.3f %.3f ", diru, dirv, dirw, eng);
    }
    fprintf(fp, "\n");
  }
  
  return;

}
