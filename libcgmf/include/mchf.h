/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file mchf.h

  \brief Monte Carlo Hauser-Feshbach class and prototype definitions

*/

#ifndef __MCHF_H__
#define __MCHF_H__

const int MAX_CASCADE    =    100 ;  /* maximum number of gamma-rays per cascade */
const int MAX_SIMULATION =      1 ;  /* number of Monte Carlo Simulations */

/****************************/
/*    Monte Carlo Index     */
/****************************/
class MCPoint{
 private:
  int    n;    // nucleus index
  int    p;    // CN state parity
  int    j;    // J-array index
  int    k;    // excitation index
  double e;    // excitation energy
  bool   l;    // true if discrete level
 public:
  MCPoint(){
    n  = 0;
    p  = 0;
    j  = 0;
    k  = 0;
    e  = 0.0;
    l  = false;
  }
  MCPoint(int n0, int p0, int j0, int k0, double e0){
    n = n0;
    p = p0;
    j = j0;
    k = k0;
    e = e0;
    l = false;
  }
  void set(int n0, int p0, int j0, int k0, double e0){
    n = n0;
    p = p0;
    j = j0;
    k = k0;
    e = e0;
  }
  void setE(double e0){
    e = e0;
  }
  void setDiscrete(){
    l = true;
  }
  int    getN(){return n;};
  int    getK(){return k;};
  int    getP(){return p;};
  int    getJ(){return j;};
  double getE(){return e;};
  bool   ifDiscrete(){return l;};
};


/****************************/
/*    Particle History      */
/****************************/
class MCHistory{
 private:
  int     id;  // particle identifier
  double  ep;  // emission energy
  MCPoint mp;  // MC index for intermediate state
 public:
  MCHistory(){
    id = 0;
    ep = 0.0;
    mp.set(0,0,0,0,0.0);
  }
  MCHistory(int id0, double ep0, MCPoint mp0){
    id = id0;
    ep = ep0;
    mp = mp0;
  }
  void set(int id0, double ep0, MCPoint mp0){
    id = id0;
    ep = ep0;
    mp = mp0;
  }
  int     getIndex()  {return id;};
  double  getEnergy() {return ep;};
  MCPoint getMCPoint(){return mp;};
};

#endif //__MCHF_H__
