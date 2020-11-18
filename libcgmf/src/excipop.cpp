/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file excipop.cpp

  \brief Calculation of emission spectra and population in exciton model

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "exciton.h"


/**********************************************************/
/*      Occupation Probability of Each p-h Configuration  */
/**********************************************************/
void  preqPOccupation(double **p, double **t1, double **t2,
                      double **tzp,double **tnp, double **tz0, double **tn0)
                             
{
  double x1,x2,x3,x4,x5,x6,y1,y2;
  bool   conv    = true;

  const double eps     = 1e-6;
  const int    maxloop = 100;
  
  for(int loop=1 ; loop<maxloop ; loop++){
    conv=true;

    for(int i=0 ; i<MAX_PREEQ ; i++){
      for(int j=0 ; j<MAX_PREEQ ; j++){

        p[0][0] = 1.0;

        x1 = x2 = x3 = x4 = x5 = x6 = y1= y2 = 0.0;

        if(i>0){
          x1 = p[i-1][j] * tzp[i-1][j] * t1[i-1][j];
          x5 = p[i-1][j] * tnp[i-1][j] * t2[i-1][j];
          if(j<MAX_PREEQ-1) y1 = tn0[i-1][j+1] * t1[i-1][j+1];
        }

        if(j>0){
          x2 = p[i][j-1] * tnp[i][j-1] * t1[i][j-1];
          x6 = p[i][j-1] * tzp[i][j-1] * t2[i][j-1];
          if(i<MAX_PREEQ-1) y2 = tz0[i+1][j-1] * t1[i+1][j-1];
        }

        if(i>1 && j<MAX_PREEQ-1) x3 = p[i-2][j+1] * tzp[i-2][j+1]* t2[i-2][j+1];
        if(j>1 && i<MAX_PREEQ-1) x4 = p[i+1][j-2] * tnp[i+1][j-2]* t2[i+1][j-2];

        double t = x1 + x2 + y1*(x3+x5) + y2*(x4+x6);
        if(t>0.0 && abs(p[i][j]/t-1.0)>eps) conv=false;
        p[i][j] = t;
      }
    }
/*
    for(int i=0 ; i<MAX_PREEQ ; i++){
      for(int j=0 ; j<MAX_PREEQ ; j++){
        cout << " " << p[i][j];
      }
      cout <<endl;
    }
    cout <<endl;
*/
    if(conv) break;
  }
  if(!conv) cout << "P calc. not converged" << endl;
}
