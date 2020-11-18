/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file ripl2levels.cpp

  \brief Read and modify nuclear structure data from RIPL

*/

#include <string>
#include <sstream>
#include <ostream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <math.h>

using namespace std;

#include "structur.h"
#include "ripl2levels.h"
#include "terminate.h"
#include "config.h"

static void riplSelectCandidate  (string, double *, int *);
static int  riplFixTransition    (int *, double *, double *, int k , Level * );
static int  riplNormalizeBranch  (int, double *, double * );
static int  riplSetNmax          (MaxLevelCtl, int, int, int);
static int  riplFixNmax          (int, Level *);
static int  riplGSonly           (unsigned int, unsigned int, Level *);

extern Level **allLevels; // all discrete levels from RIPL database
extern int *zaidIndex;
extern int *numberLevels;
static int NUM_NUCLEI = 0;
const int MAX_NUCLEI = 3500;

/*!
 
 \brief Read complete database of CGMF discrete levels, based on RIPL3 file.
 
  HEADER
  ======

 This routine reads the CGMF-ready discrete level file, which has been processed
 for use in CGMF already. It does not read the RIPL3 database directly. The 
 advantage is that the entire file is 4.8 MB instead of ~60 MB for RIPL3. This 
 file size can be reduced even more if considering only fission fragments, but
 this would be a problem for CGM (non-fission events). The disadvantage is that
 not all information contained in RIPL3 is kept in the CGMF data file.
 
 */
void readDiscreteLevelData () {
	
	ostringstream os;
	ifstream		  fp;
	string			  str,file,d;
	
	int zaid;
	
	int nlevel = 0;
	
	MaxLevelCtl ml; // to define, used to be in the argument list
	ml=reassign;

	ZAnumber za;
	
	int 	*nf;
	double 	*pe, *ic;
	
	allLevels    = new Level * [MAX_NUCLEI];
	zaidIndex    = new int [MAX_NUCLEI];
	numberLevels = new int [MAX_NUCLEI];
	
	for (int i=0; i<MAX_NUCLEI; i++) {
		zaidIndex[i]=0;
		numberLevels[i]=0;
	}

	str = DISCRETEDATA;
	
	// try current directory first
	fp.open(&str[0]);
	
	// then system data area
	if (!fp) {
		str = datadir + DISCRETEDATA;
		fp.open(&str[0]);
	}
	
	if(!fp){
		cgmTerminateCode("CGMF discrete level database not found");
		return;
	}
	
	int a=0, z=0, nol=0, nc=0, nmax=0;
	int c=0;

	while(getline(fp,str)){ // loop over data file
		
		if (str[0]=='#') continue;
		
		d=str.substr(6,6); zaid=atoi(&d[0]); // ZAID
		d=str.substr(12,5); nol=atoi(&d[0]); // number of levels

		if(nol >= MAX_LEVELS){
			cerr << "discrete levels for Z " << za.getZ() << " - A " << za.getA()
			<< " too many (" << nol << ")";
			fp.close();
			cgmTerminateCode(os.str());
		}
		
		zaidIndex[c] = zaid;
		numberLevels[c] = nol;
		allLevels[c] = new Level [nol];
		
		/*** for all discrete levels */
		double elev=0.0, s=0.0, thlf=0.0;
		int    p=0;
		int    ng=0;
		
		for(int i=0 ; i<nol ; i++){ // loop over discrete levels
			
			getline(fp,str);
			d=str.substr( 0,10);  elev = atof(&d[0]); // energy (MeV)
			d=str.substr(10, 5);  s    = atof(&d[0]); // spin (hbar)
			d=str.substr(15, 3);  p    = atoi(&d[0]); // parity
			d=str.substr(18, 13); thlf = atof(&d[0]); // half-life (sec)
			d=str.substr(31, 3);  ng   = atoi(&d[0]); // number of gamma rays de-exciting the level

			nf = new int [max(ng,MAX_GAMMA_BRANCH)];
			pe = new double [max(ng,MAX_GAMMA_BRANCH)];
			ic = new double [max(ng,MAX_GAMMA_BRANCH)];

			for(int j=0 ; j<max(ng,MAX_GAMMA_BRANCH) ; j++){
				nf[j] = 0;
				pe[j] = 1.0;
				ic[j] = 0.0;
			}
			
			for(int j=0 ; j<ng ; j++){
				getline(fp,str);
				d=str.substr( 0, 3); nf[j] = atoi(&d[0]);   // serial number of the final state
				d=str.substr( 3, 9); pe[j] = atof(&d[0]);   // prob. of decay by electromagnetic transition
				d=str.substr(12,10); ic[j] = atof(&d[0]);
			}


				// if( (s < 0.0) && (i == 0) ) {
				// 	riplSelectCandidate(str, &s, &p); // ground-state spin unknown
				// }
				
			/* moved this above next block, otherwise spin,parity of last level is set to zero! */
			allLevels[c][i].energy   = elev;
			allLevels[c][i].spin     = s;
			allLevels[c][i].halflife = thlf;
			allLevels[c][i].parity   = p;
				
			if(ng==0){ // no gamma assigned
				allLevels[c][i].fstate[0] = 0;
				allLevels[c][i].branch[0] = 1.0;
				allLevels[c][i].gratio[0] = 1.0;
				ng=riplFixTransition(nf,pe,ic,i,allLevels[c]);
			}
			ng = riplNormalizeBranch(ng,pe,ic);
				
			allLevels[c][i].ngamma   = ng;
				
			for(int j=0 ; j<min(ng,MAX_GAMMA_BRANCH) ; j++){
				allLevels[c][i].fstate[j] = nf[j];
				allLevels[c][i].branch[j] = pe[j];
				allLevels[c][i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
			}

			delete [] nf;
			delete [] pe;
			delete [] ic;

		} // end loop over discrete levels

		c++;

	}
	
	fp.close();
	
	NUM_NUCLEI = c;
	
	return;
	
}

/*! \brief Retrieve CGMF discrete levels for a specific nucleus
 
 Retrieve CGMF discrete levels for a specific nucleus, and the
 ones that follow from its decay up to nemit neutrons, once
 the entire database is read using cgmfReadDiscreteLevelData.

 If a nucleus is not found, default ground-state spin and parity 
 are set by calling the riplGSonly routine.
 
 */
void getDiscreteLevels (Nucleus *ncl, int nemit) {
		
	int id;
	int zaid;
	bool found;

	for (int k=0; k<=nemit; k++) {
		zaid = 1000*ncl[k].za.getZ()+ncl[k].za.getA();
		found=false;
		id=-1;
		while (!found && id<MAX_NUCLEI) {
			id++;
			if (zaidIndex[id]==zaid) found=true;
		}
		if (!found) {
			ncl[k].ndisc = riplGSonly(ncl[k].za.getZ(), ncl[k].za.getA(), ncl[k].lev);
		} else {
			ncl[k].ndisc = numberLevels[id];
			for (int i=0; i<ncl[k].ndisc; i++) {
				ncl[k].lev[i] = allLevels[id][i];
			}
		}
	}

	return;

}


/*! \brief Retrieve CGMF discrete levels for a specific nucleus.

Retrieve the discrete levels for a specific nucleus from a CGMF data file, 
which has processed the RIPL3 database first.

*/
void getDiscreteLevels (Nucleus *ncl) {
	cout << "Here?\n";
	getDiscreteLevels (ncl, 0);
	return;
}

/*! \brief Retrieve RIPL discrete levels for a specific nucleus
 
 Retrieve RIPL discrete levels for a specific nucleus, and the
 ones that follow from its decay up to nemit neutrons, once
 the entire database is read using riplReadDiscreteLevelsData
 
 */
void getRiplDiscreteLevels (Nucleus *ncl, int nemit) {
	
	MaxLevelCtl ml = reassign;
	
	int zaid = 1000*ncl[0].za.getZ()+ncl[0].za.getA();

	bool found =false;
	int id=-1;
	while (!found && id<MAX_NUCLEI-1) {
		id++;
		if (zaidIndex[id]==zaid) found=true;
	}
	
	for (int k=0; k<=nemit; k++) {
		found = false;
		if (zaidIndex[id-k]==zaid-k) found = true;
		ncl[k].ndisc = numberLevels[id-k];
		for (int i=0; i<ncl[k].ndisc; i++) {
			ncl[k].lev[i] = allLevels[id-k][i];
		}
		fixLevels (&ncl[k], found, ml);
	}
	
	return;
	
}

/*!
 \brief Retrieve the discrete levels for a specific nucleus.
 \param *ncl: pointer to a Nucleus
 */
void getRiplDiscreteLevels (Nucleus *ncl) {
	getRiplDiscreteLevels(ncl, 0);
	return;
}

/*!
 
 \brief Fix levels that are missing spin or/and parity assignements.
 
 Discrete levels are read from the RIPL database. Some levels lack proper
 spin and/or parity assignements, and need to be fixed before CGMF can use
 them.
 
 \param *ncl: pointer to a Nucleus
 \param found: boolean stating that the nucleus has been found in the data file
 \param ml: flag indicating how to treat levels with missing spin and/or parity
 
 \todo I don't think I need this routine anymore once I finalize my new level file for CGMF.
 
 */
void fixLevels (Nucleus *ncl, bool found, MaxLevelCtl ml) {
	
	ostringstream os;
	
	if(!found){
		if(ml != reassign){
			os.str("");
			os << "discrete level data for Z " << ncl[0].za.getZ() << " - A " <<
				ncl[0].za.getA() << " not found";
			cgmTerminateCode(os.str());
		}
		else{
			ncl[0].ndisc = riplGSonly(ncl[0].za.getZ(), ncl[0].za.getA(), ncl[0].lev);
		}
	}
	else{
		if( (ncl[0].lev[0].spin < 0.0) || ncl[0].lev[0].parity == 0 ){
			if(ml == reassign)
				ncl[0].ndisc = riplGSonly(ncl[0].za.getZ(),ncl[0].za.getA(), ncl[0].lev);
			else{
				os.str("");
				os << "discrete level data for Z " << ncl[0].za.getZ() << " - A " <<
					ncl[0].za.getA() << " found but spin/parity not assigned";
				cgmTerminateCode(os.str());
			}
		}
	}
	
	if(ml == extended){
		ncl[0].ndisc = riplFixNmax(ncl[0].ndisc,ncl[0].lev);
		if(ncl[0].ndisc == 0){
			os.str("");
			os << "discrete level data for Z " << ncl[0].za.getZ() << " - A " <<
				ncl[0].za.getA() << " not complete";
			cgmTerminateCode(os.str());
		}
	}
	
}


/*!
 
 \brief Read entire file of discrete levels from RIPL database.
 
 \todo Rewrite this routine to read new CGMF level file
 
 \attention THIS IS THE ONE I AM USING!
 
 \todo Clean all routines not used anymore!
 
 */
void riplReadDiscreteLevelData () {

	ostringstream os;
	ifstream		  fp;
	string			  str,file,d;

	int nlevel = 0;
	
	MaxLevelCtl ml; // to define, used to be in the argument list
	ml=reassign;
	
	ZAnumber za;
	
	int           nf[MAX_GAMMA_BRANCH];
	double        pe[MAX_GAMMA_BRANCH], ic[MAX_GAMMA_BRANCH];
	
	allLevels    = new Level * [MAX_NUCLEI];
	zaidIndex    = new int [MAX_NUCLEI];
	numberLevels = new int [MAX_NUCLEI];

	str = RIPLDATA;

	// try current directory first
	fp.open(&str[0]);
	
	// then system data area
	if (!fp) {
		str = datadir + RIPLDATA;
		fp.open(&str[0]);
	}
	
	if(!fp){
		cgmTerminateCode("discrete level database not found");
		return;
	}
	
	int a=0, z=0, nol=0, nc=0, nmax=0;
	int c=0;
	
	while(getline(fp,str)){ // loop over data file

		for(int i=0 ; i<MAX_GAMMA_BRANCH ; i++) {
			nf[i] = 0;
			pe[i] = 1.0;
			ic[i] = 0.0;
		}	
		
		/*** search for Z and A entry in the file */
		d=str.substr( 5, 5);  a    = atoi(&d[0]); // mass number
		d=str.substr(10, 5);  z    = atoi(&d[0]); // charge number
		d=str.substr(15, 5);  nol  = atoi(&d[0]); // number of levels in decay scheme
		d=str.substr(20, 5);//nog  = atoi(&d[0]); // number of gamma rays in decay scheme
		d=str.substr(25, 5);  nmax = atoi(&d[0]); // max. number of levels up to which the level scheme is complete
		d=str.substr(30, 5);  nc   = atoi(&d[0]); // level number up to which the spins and parities are unique
		ZAnumber za(z,a);
		
		zaidIndex[c] = 1000*z+a;

		if(nmax >= MAX_LEVELS){
			cerr << "discrete levels for Z " << za.getZ() << " - A " << za.getA()
				<< " too many (" << nol << ")";
			fp.close();
			cgmTerminateCode(os.str());
		}
		
		/*** set Nmax by the value of MaxLevelCtl */
		nlevel = riplSetNmax(ml, nol, nmax, nc);
		numberLevels[c] = nlevel;
		allLevels[c] = new Level [nlevel]; // changed from nol to nlevel

		/*** for all discrete levels */
		double elev=0.0, s=0.0, thlf=0.0;
		int    p=0;
		int    ng=0;

		for(int i=0 ; i<nol ; i++){ // loop over discrete levels

			getline(fp,str);
			d=str.substr( 4,10);  elev = atof(&d[0]); // energy (MeV)
			d=str.substr(15, 5);  s    = atof(&d[0]); // spin (hbar)
			d=str.substr(20, 3);  p    = atoi(&d[0]); // parity
			d=str.substr(25, 9);  thlf = atof(&d[0]); // half-life (sec)
			d=str.substr(34, 3);  ng   = atoi(&d[0]); // number of gamma rays de-exciting the level
			
			for(int j=0 ; j<ng ; j++){
				getline(fp,str);
				d=str.substr(39, 4); nf[j] = atoi(&d[0]) -1; // serial number of the final state
				d=str.substr(66,10); pe[j] = atof(&d[0]);    // prob. of decay by electromagnetic transition
				d=str.substr(77,10); ic[j] = atof(&d[0]);
			}
			
			if (i<nlevel) { // store only data up to nlevel
				
				if( (s < 0.0) && (i == 0) ) {
					riplSelectCandidate(str, &s, &p); // ground-state spin unknown
				}

				/* moved this above next block, otherwise spin,parity of last level is set to zero! */
				allLevels[c][i].energy   = elev;
				allLevels[c][i].spin     = s;
				allLevels[c][i].halflife = thlf;
				allLevels[c][i].parity   = p;

				if(ng==0){ // no gamma assigned
					allLevels[c][i].fstate[0] = 0;
					allLevels[c][i].branch[0] = 1.0;
					allLevels[c][i].gratio[0] = 1.0;
					ng=riplFixTransition(nf,pe,ic,i,allLevels[c]);
				}
				ng = riplNormalizeBranch(ng,pe,ic);
				
				allLevels[c][i].ngamma   = ng;

				for(int j=0 ; j<min(ng,MAX_GAMMA_BRANCH) ; j++){
					allLevels[c][i].fstate[j] = nf[j];
					allLevels[c][i].branch[j] = pe[j];
					allLevels[c][i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
				}
				
			}
			
		} // end loop over discrete levels
		
		c++;

	}

	fp.close();

	NUM_NUCLEI = c;

	return;
	
}



/**********************************************************/
/*      Read Discrete Level Data from RIPL Database       */
/**********************************************************/
int riplReadDiscreteLevels(ZAnumber *za, Level *lev, MaxLevelCtl ml)
{
  ostringstream os;
  ifstream      fp;
  string        str,file,d;
  int           nlevel = 0;
  int           nf[MAX_GAMMA_BRANCH];
  double        pe[MAX_GAMMA_BRANCH],ic[MAX_GAMMA_BRANCH];
  
  os << setw(4) << setfill('0') << za->getZ();

  file = os.str();
  file[0] = 'z';
  str = file + ".dat";
  
  /*** try current directry first */
  fp.open(&str[0]);
  
  /*** then system data area */
  if(!fp){
    str = datadir + LEVELDIRECTORY + file + ".dat";
    fp.open(&str[0]);
  }
  
  if(!fp){
    cgmTerminateCode("discrete level data file not found");
    nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    return(nlevel);
  }
	
  for(int i=0 ; i<MAX_LEVELS ; i++) lev[i].energy = lev[i].spin = 0.0;
  
  bool found=false;
  int  a=0, z=0, nol=0, nc=0, nmax=0;
  while(getline(fp,str)){
    
    /*** search for Z and A entry in the file */
    d=str.substr( 5, 5);  a    = atoi(&d[0]);
    d=str.substr(10, 5);  z    = atoi(&d[0]);
    d=str.substr(15, 5);  nol  = atoi(&d[0]);
    d=str.substr(20, 5);//nog  = atoi(&d[0]);
    d=str.substr(25, 5);  nmax = atoi(&d[0]);
    d=str.substr(30, 5);  nc   = atoi(&d[0]);
    ZAnumber za1(z,a);
    
    if((za1.getZ()==za->getZ()) && (za1.getA()==za->getA())) found = true;
    
    if(nmax >= MAX_LEVELS){
      os.str("");
      os << "discrete levels for Z " << za1.getZ() << " - A " << za1.getA() << " too many (" << nol << ")";
      fp.close();
      cgmTerminateCode(os.str());
    }
    
    /*** set Nmax by the value of MaxLevelCtl */
    nlevel = riplSetNmax(ml, nol, nmax, nc);
        
    /*** for all discrete levels */
    double elev=0.0, s=0.0, thlf=0.0;
    int    p=0, ng=0;
    for(int i=0 ; i<nol ; i++){
      getline(fp,str);
      d=str.substr( 4,10);  elev = atof(&d[0]);
      d=str.substr(15, 5);  s    = atof(&d[0]);
      d=str.substr(20, 3);  p    = atoi(&d[0]);
      d=str.substr(25, 9);  thlf = atof(&d[0]);
      d=str.substr(34, 3);  ng   = atoi(&d[0]);
            
      /*** if no gamma is assigned */
      if(ng==0){
        nf[0] = 0;
        pe[0] = 1.0;
        ic[0] = 0.0;
      }
      /*** for gamma-ray branches */
      else{
        for(int j=0 ; j<ng ; j++){
          getline(fp,str);
          d=str.substr(39, 4);  nf[j] = atoi(&d[0]) -1;
          d=str.substr(66,10);  pe[j] = atof(&d[0]);
          d=str.substr(77,10);  ic[j] = atof(&d[0]);
        }
      }
      
      /*** if found, copy data into lev array */
      if(found && i<nlevel){
        
        if( (s < 0.0) && (i == 0) ) riplSelectCandidate(str, &s, &p);
        if(ng==0){
          ng=riplFixTransition(nf,pe,ic,i,lev);
        }
        ng = riplNormalizeBranch(ng,pe,ic);
        
        lev[i].energy   = elev;
        lev[i].spin     = s;
        lev[i].halflife = thlf;
        lev[i].parity   = p;
        lev[i].ngamma   = ng;
        
        for(int j=0 ; j<ng ; j++){
          lev[i].fstate[j] = nf[j];
          lev[i].branch[j] = pe[j];
          lev[i].gratio[j] = 1.0/(1.0 + ic[j]); // Ig/(Ig + Ie) = 1/(1+ICC)
        }
      }
    }
    if(found) break;
  }
  fp.close();
    
  if(!found){
    if(ml != reassign){
      os.str("");
      os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not found";
      cgmTerminateCode(os.str());
    }
    else{
      nlevel = riplGSonly(za->getZ(), za->getA(), lev);
    }
  }
  else{
    if( (lev[0].spin < 0.0) || lev[0].parity == 0 ){
      if(ml == reassign) nlevel = riplGSonly(za->getZ(), za->getA(), lev);
      else{
        os.str("");
        os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " found but spin/parity not assigned";
        cgmTerminateCode(os.str());
      }
    }
  }
  
  if(ml == extended){
    nlevel = riplFixNmax(nlevel,lev);
    if(nlevel == 0){
      os.str("");
      os << "discrete level data for Z " << za->getZ() << " - A " << za->getA() << " not complete";
      cgmTerminateCode(os.str());
    }
  }
    
  return(nlevel);
}



/**********************************************************/
/*      Get Spin and Parity from Candidates               */
/**********************************************************/
void riplSelectCandidate(string str, double *s, int *p)
{
  string csrc = str.substr(46,18);
  char   cdst[] = "                  ";
	
  int i0 = csrc.find("(");
  int i1 = csrc.find_first_of(",)");
	
  if(i0 == (signed int)string::npos) return;

  for(int i=i0+1 ; i<i1 ; i++) cdst[i-i0-1] = csrc[i];
  cdst[i1-i0-1] = '\0';
	
  *s = (strchr(cdst,'/') != NULL) ?  0.5 : 1.0;
  *p = (strchr(cdst,'-') != NULL) ? -1   : 1;
	
  i1 = strlen(cdst);
  for(int i=0 ; i<i1 ; i++){
    if( (cdst[i] == '+') || (cdst[i] == '-') || (cdst[i] == '/') ){
      cdst[i] = '\0';
      break;
    }
  }
  *s *= atof(cdst);

}

/**********************************************************/
/*      Fix Transition If No Decay                        */
/**********************************************************/
int riplFixTransition(int nf[], double pe[], double ic[], int k, Level *lev)
{
  int n = 0;
	
  /*** if this is the first excited state */
  if(k == 1){
    nf[0] = 0;
    pe[0] = 1.0;
    ic[0] = 0.0;
    n = 1;
    return(n);
  }
  
  /*** first, look for states to which decay is possible by E1 */
  int    p = lev[k].parity;
  double s = lev[k].spin;
	
  int dm = 100;
  for(int i=0 ; i<k ; i++){
    if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
    int dj = (int)fabs(lev[i].spin - s );
    if(dj < dm) dm = dj;
    /*** E1, different parity and |J0 - J1| <= 1.0 */
		if( (lev[i].parity != p) && (dj <= 1) ){ /*! \bug could be wrong if J0=J1 */
      nf[n] = i;
      pe[n] = 1.0;
      ic[n] = 0.0;
      n++;
      if(n == MAX_GAMMA_BRANCH-1) break;
    }
  }
	
  /*** if no transition still,
   find the closest spin state with the lowest energy */
	/*! \bug Why choose the state with lowest energy? Wouldn't it be better to choose the one closest in energy instead? */
  if(n == 0){
    for(int i=0 ; i<k ; i++){
      if( (s == 0.0) && (lev[i].spin == 0.0) ) continue;
      int dj = (int)fabs(lev[i].spin - s);
      if(dj == dm){
        nf[0] = i;
        pe[0] = 1.0;
        ic[0] = 0.0;
        n = 1;
        break;
      }
    }
  }
  
  return (n) ;
}

/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplNormalizeBranch(int ng, double pe[], double ic[])
{
  double s=0.0;
  
  for(int i=0 ; i<ng ; i++) s += pe[i];
  
  if(s == 0.0){
    if(ng > 1){
      s = 1.0/(double)ng;
      for(int i=0 ; i<ng ; i++){
        pe[i] = s;
        ic[i] = 0.0;
      }
    }
    else{
      pe[0] = 1.0;  // decay 100% to the first line in the branches
      ic[0] = 0.0;
      ng = 1;
    }
  }else{
    s = 1.0/s;
    for(int i=0 ; i<ng ; i++) pe[i] *= s;
  }
  
  return(ng);
}


/**********************************************************/
/*      Determine Nmax from Input                         */
/**********************************************************/
int riplSetNmax(MaxLevelCtl ml, int nol, int nmax, int nc)
{
  int nlevel = 0;
  
  switch(ml){
    case normal:   // Nc is used for the highest level
      nlevel = nc;
      break;
    case extended: // Nmax is used, but up to unknown spin/parity state
    case reassign: // Nmax is used, and unknown level spin/parity will be re-assigned later
      nlevel = nmax;
      break;
    case all:      // include all levels given
      nlevel = nol;
      break;
    default:
      nlevel = nc;
      break;
  }
	
	if (nlevel>MAX_LEVELS) nlevel=MAX_LEVELS; // cut at MAX_LEVELS
	
  return(nlevel);
}


/**********************************************************/
/*      Renormalize Transition Probability                */
/**********************************************************/
int riplFixNmax(int nl, Level *lev)
{
  for(int i=0 ; i<nl ; i++){
    /*** spin or parity not assigned */
    if( (lev[i].spin < 0.0) || (lev[i].parity == 0) ){
      nl = i;
      break;
    }
  }
  
  return(nl);
}


/**********************************************************/
/*      Assume Ground State Spin                          */
/**********************************************************/
int riplGSonly(unsigned int z, unsigned int a, Level *lev)
{
  unsigned int n = a -z;
  
  /*** even-even */
  if(((z%2) == 0) && ((n%2) == 0) ){
    lev[0].spin   = 0.0;
    lev[0].parity = 1;
  }
  /*** odd-odd, 1+ assumed */
  else if(((z%2) != 0) && ((n%2) != 0) ){
    lev[0].spin   = 1.0;
    lev[0].parity = 1;
  }
  /*** even-odd, 1/2+ assumed */
  else{
    lev[0].spin   = 0.5;
    lev[0].parity = 1;
  }
  
  return(1);
}



void riplDiscreteLevelsCleanup()
{

	for (int i=0; i<NUM_NUCLEI; i++) {
		delete [] allLevels[i];
	}

	delete [] allLevels;
	delete [] zaidIndex;
	delete [] numberLevels;
  
}
