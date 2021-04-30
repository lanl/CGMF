/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file excilambda.cpp

  \brief Lambda functions for lifetime calculation

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "physics.h"
#include "exciton.h"

static inline double  preqLpz0  (Preeq *, Exconf *);
static inline double  preqLpn0  (Preeq *, Exconf *);
static inline double  densinteg (double, Preeq *, Exconf *, Exconf *, Exconf *);
static inline double  epauli    (int, int, int, int, Exconf *, Preeq *);

/**********************************************************/
/*      Internal Transition Rate for Proton               */
/*      All rates calculated here are in h-bar unit.      */
/*      Dividing by 6.58E-22 gives the absolute rate [/s] */
/**********************************************************/
double preqLambdaPlusZ(Preeq *q, Exconf *e)
{
 if( (e->zp + e->zh + e->np + e->nh)==1 ) return(preqLpz0(q,e));

 double ep[4];
 Exconf xp[4],xr[4],xw[4];

 ep[0] = epauli(-1, 0, 0, 0,e,q);
 ep[1] = epauli( 0,-1, 0, 0,e,q);
 ep[2] = epauli( 0, 0,-1, 0,e,q);
 ep[3] = epauli( 0, 0, 0,-1,e,q);

 xp[0].setExconf(      1,    0,    0,    0);
 xr[0].setExconf(e->zp-1,e->zh,e->np,e->nh);
 xw[0].setExconf(      2,    1,    0,    0);

 xp[1].setExconf(    0,      1,    0,    0);
 xr[1].setExconf(e->zp,e->zh-1,e->np,e->nh);
 xw[1].setExconf(    1,      2,    0,    0);

 xp[2].setExconf(    0,    0,      1,    0);
 xr[2].setExconf(e->zp,e->zh,e->np-1,e->nh);
 xw[2].setExconf(    1,    1,      1,    0);

 xp[3].setExconf(    0,    0,    0,      1);
 xr[3].setExconf(e->zp,e->zh,e->np,e->nh-1);
 xw[3].setExconf(    1,    1,    0,      1);

 double c1 = PI2 * q->m2zz;
 double c2 = PI2 * q->m2nz;

 const int ndiv = 100;
 double lambda  = 0.0;

 for(int k=0 ; k<4 ; k++){
   double e1 = epauli(1,1,0,0,e,q) - ep[k];
   double e2 = q->ex_total         - ep[k];
   double de = (e2 - e1)/(double)ndiv;
   double f0 = de/3.0;
   double f1 = 2.0*f0;
   double f2 = 4.0*f0;

   lambda += ( densinteg(e1,q,&xw[k],&xp[k],&xr[k])
              +densinteg(e2,q,&xw[k],&xp[k],&xr[k]) )*f0;

   for(int i=1 ; i<ndiv ; i++){
     lambda += ((i%2==0) ? f1 : f2) * ((k<2) ? c1 : c2)
              *densinteg(e1+i*de,q,&xw[k],&xp[k],&xr[k]);
   }
 }

 return( lambda/q->omega_total );
}


/**********************************************************/
/*      Internal Transition Rate for Neutron              */
/**********************************************************/
double preqLambdaPlusN(Preeq *q, Exconf *e)
{
 if( (e->zp + e->zh + e->np + e->nh)==1 ) return(preqLpn0(q,e));

 double ep[4];
 Exconf xp[4],xr[4],xw[4];

 ep[0] = epauli( 0, 0,-1, 0,e,q);
 ep[1] = epauli( 0, 0, 0,-1,e,q);
 ep[2] = epauli(-1, 0, 0, 0,e,q);
 ep[3] = epauli( 0,-1, 0, 0,e,q);

 xp[0].setExconf(    0,    0,      1,    0);
 xr[0].setExconf(e->zp,e->zh,e->np-1,e->nh);
 xw[0].setExconf(    0,    0,      2,    1);

 xp[1].setExconf(    0,    0,    0,      1);
 xr[1].setExconf(e->zp,e->zh,e->np,e->nh-1);
 xw[1].setExconf(    0,    0,    1,      2);

 xp[2].setExconf(      1,    0,    0,    0);
 xr[2].setExconf(e->zp-1,e->zh,e->np,e->nh);
 xw[2].setExconf(      1,    0,    1,    1);

 xp[3].setExconf(    0,      1,    0,    0);
 xr[3].setExconf(e->zp,e->zh-1,e->np,e->nh);
 xw[3].setExconf(    0,      1,    1,    1);

 double c1 = PI2 * q->m2nn;
 double c2 = PI2 * q->m2zn;

 const int ndiv = 100;
 double lambda  = 0.0;

 for(int k=0 ; k<4 ; k++){
   double e1 = epauli(0,0,1,1,e,q) - ep[k];
   double e2 = q->ex_total         - ep[k];
   double de = (e2 - e1)/(double)ndiv;
   double f0 = de/3.0;
   double f1 = 2.0*f0;
   double f2 = 4.0*f0;

   lambda += ( densinteg(e1,q,&xw[k],&xp[k],&xr[k])
              +densinteg(e2,q,&xw[k],&xp[k],&xr[k]) )*f0;

   for(int i=1 ; i<ndiv ; i++){
     lambda += ((i%2==0) ? f1 : f2) * ((k<2) ? c1 : c2)
              *densinteg(e1+i*de,q,&xw[k],&xp[k],&xr[k]);
   }
 }

 return( lambda/q->omega_total);
}


/**********************************************************/
/*      Internal Transision (special case for n=1)        */
/**********************************************************/
double preqLpz0(Preeq *q, Exconf *e)
{
 double lambda = 0.0;

 double e2 = q->ex_total - epauli(1,1,0,0,e,q);

 if(e2>0.0){
   double c  = PI2/4.0 * q->spd[0].gz*q->spd[0].gz;
   double c1 =     (e->zp+e->zh) * q->spd[0].gz * q->m2zz;
   double c2 = 2.0*(e->np+e->nh) * q->spd[0].gn * q->m2zn;
   double c3 = preqFiniteWell(3,1,q->ex_total,q->spd[0].well_depth);

   lambda = c * e2*e2 * (c1+c2) * c3;
 }
 return(lambda);
}


double preqLpn0(Preeq *q, Exconf *e)
{
 double lambda = 0.0;

 double e2 = q->ex_total - epauli(0,0,1,1,e,q);

 if(e2>0.0){
   double c  = PI2/4.0 * q->spd[0].gn*q->spd[0].gn;
   double c1 =     (e->np+e->nh) * q->spd[0].gn * q->m2nn;
   double c2 = 2.0*(e->zp+e->zh) * q->spd[0].gz * q->m2nz;
   double c3 = preqFiniteWell(3,1,q->ex_total,q->spd[0].well_depth);

   lambda = c * e2*e2 * (c1+c2) * c3;
 }
 return(lambda);
}


/**********************************************************/
/*      Conversion Rate for Proton                        */
/**********************************************************/
double preqLambdaZeroZ(Preeq *q, Exconf *e)
{
 double ep = epauli(-1,-1,0,0,e,q);
 double e1 = epauli( 0, 0,0,0,e,q) - ep;
 double e2 = q->ex_total           - ep;

 double c1 = PI2 * q->m2zn;
 Exconf xp,xr,xw;

 xp.setExconf(      1,      1,    0,    0);
 xr.setExconf(e->zp-1,e->zh-1,e->np,e->nh);
 xw.setExconf(      0,      0,    1,    1);

 const int ndiv = 100;
 double lambda  = 0.0;

 double de = (e2 - e1)/(double)ndiv;
 double f0 = de/3.0;
 double f1 = 2.0*f0;
 double f2 = 4.0*f0;

 lambda += ( densinteg(e1,q,&xw,&xp,&xr)
            +densinteg(e2,q,&xw,&xp,&xr) )*f0;

 for(int i=1 ; i<ndiv ; i++){
   lambda += ((i%2==0) ? f1 : f2) * c1
            * densinteg(e1+i*de,q,&xw,&xp,&xr);
 }

 return( lambda/q->omega_total );
}


/**********************************************************/
/*      Conversion Rate for Neutron                       */
/**********************************************************/
double preqLambdaZeroN(Preeq *q, Exconf *e)
{
 double ep = epauli(0,0,-1,-1,e,q);
 double e1 = epauli(0,0, 0, 0,e,q) - ep;
 double e2 = q->ex_total           - ep;

 double c1 = PI2 * q->m2nz;
 Exconf xp,xr,xw;

 xp.setExconf(    0,    0,      1,      1);
 xr.setExconf(e->zp,e->zh,e->np-1,e->nh-1);
 xw.setExconf(    1,    1,      0,      0);

 const int ndiv = 100;
 double lambda  = 0.0;

 double de = (e2 - e1)/(double)ndiv;
 double f0 = de/3.0;
 double f1 = 2.0*f0;
 double f2 = 4.0*f0;

 lambda += ( densinteg(e1,q,&xw,&xp,&xr)
            +densinteg(e2,q,&xw,&xp,&xr) )*f0;

 for(int i=1 ; i<ndiv ; i++){
   lambda += ((i%2==0) ? f1 : f2) * c1
            * densinteg(e1+i*de,q,&xw,&xp,&xr);
 }

 return( lambda/q->omega_total );
}


/**********************************************************/
/*      Integrant 3-State Densities                       */
/**********************************************************/
double densinteg(double x, Preeq *q, Exconf *ew, Exconf *ep, Exconf *er)
{
 double z = preqStateDensity(x,&q->spd[0],ew)
          * preqStateDensity(x,&q->spd[0],ep)
          * preqStateDensity(q->ex_total-x,&q->spd[0],er);
 return(z);
}


/**********************************************************/
/*      Integrant Energy Shift                            */
/**********************************************************/
double epauli(int dzp, int dzh, int dnp, int dnh, Exconf *e, Preeq *q)
{
 /*** State density for the composite system appears in the first element */
 double z = preqPauliCorrection(e->zp + dzp, e->zh + dzh,
                                e->np + dnp, e->nh + dnh,
                                q->spd[0].gz,q->spd[0].gn);
 return(z);
}
