/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file asympt.cpp

  \brief Asymptotic Coulomb functions in free space

*/

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

#include "physics.h"
#include "optical.h"
#include "maths.h"

static const double MESH_RHO   =     0.1;
static const double CRIT_CNVRG = 1.0e-06;
static const double BETA       =    0.25;

static void omGcase1 (double, double  , double *, double *);
static void omGcase2 (double, double *, double *);
static void omGcase3 (double, double  , double  , double *, double *);
static void omGcase4 (double, double  , double *, double *);
static void omNewmark(int   , double  , double  , double  , double *, double *);


/***********************************************************/
/*      C.E.Froeberg                                       */
/*      Rev. Mod. Phys., 27, 399(1955)                     */
/***********************************************************/
void omAsymptotic(double rm, double yeta, double sigma0, double *g0, double *g1)
{
   double ra,mesh,dr,p0,p1;
   int n;

   if(yeta>=15.0){
     ra=16.0;
     if(ra<rm){
       if(rm<=5.0*yeta/3.0-5.0){
         ra=rm;
         omGcase1(ra,yeta,&p0,&p1);
       }
       else{
         ra=2.0*yeta;
         omGcase2(yeta,&p0,&p1);
       }
     }
     else
       omGcase1(ra,yeta,&p0,&p1);
   }
   else if(yeta>=5.0){
     ra=2.0*yeta;
     if(ra<rm){
       if((ra=2.5*yeta+7.5)<rm) ra=rm;
       omGcase4(ra,yeta,&p0,&p1);
     }
     else
       omGcase2(yeta,&p0,&p1);
   }
   else{
     if((ra=2.0*yeta+10.0)<rm) ra=rm;
     omGcase3(ra,yeta,sigma0,&p0,&p1);
   }

   if(p0!=0.0 && p1!=0.0){
     if((dr=ra-rm)!=0.0){
       n=abs((int)( dr/MESH_RHO )+1);
       mesh=dr/(double)n;
       omNewmark(n,mesh,yeta,ra,&p0,&p1);
     }
   }
   *g0=p0;
   *g1=p1;
}

/***************************************/
/*     Riccati - 1                     */
/***************************************/
void omGcase1(double ra, double y, double *p0, double *p1)
{
   double t,t1,t2,t3,y1,g0,g1,g2,g3,g4,g5,g6,g7,q;

   y1= 1.0/(2.0*y);
   t = ra*y1;
   t1= 1.0-t;
   t2= t*t1*t1*t1;
   t3= sqrt(t2);

   g0=  sqrt(t*t1)+asin(sqrt(t))-PI_2;
   g1=  log(t/t1)/4.0;
   g2= -(8.0*t*t-12.0*t+9.0)/(48.0*t3);
   g3=  (8.0*t-3.0)/(64.0*t2);
   g4=  (t*(t*(t*(t*(t*(2048.*t-9216.)+16128.)-13440.)-12240.)+7560.)-1890.)
       /(   92160.*t2   *t3);
   g5=  3*(t*(t*(1024.*t-448.)+208.)-39.)/(8192.*t2*t2);
   g6= -(t*(t*(t*(t*(t*(t*(t*(t*(t*(262144.*t-1966080.)+6389760.)-11714560.)+
         13178880.)-9225216.)+13520640.)-3588480.)+2487240.)-873180.)+130977.)
       /(10321920.*t2*t2*t3);
   g7=  (t*(t*(t*(t*(1105920.*t-55296.)+314624.)-159552.)+45576.)-5697.)
       /(  393216.*t2*t2*t2);

   q = -2.0*y*g0+g1-y1*(g2-y1*(g3-y1*(g4-y1*(g5-y1*(g6-g7*y1)))));

/* to avoid overflow at very low energies */
   if(q>100.0){
      *p0 = *p1 = 0.0;
     return;
   }

  *p0=  exp(q);

   g0=  sqrt(t1/t);
   g1=  1.0/(4.0*t*t1);
   g2= -(8.0*t-3.0)/(32.0*sqrt(pow(t,3.0)*pow(t1,5.0)));
   g3=  3.0*(8.0*t*t-4.0*t+1.0)/(64.0*t*t*pow(t1,4.0));
   g4= -(t*(t*(1536.*t-704.)+336.)-63.)/(2048.*sqrt(pow(t,5.0)*pow(t1,11.0)));
   g5=  3.0*(t*(t*(t*(2560.*t-832.)+728.)-260.)+39.)/(4096.*t*t*t*pow(t1,7.0));
   g6=  (t*(t*(t*(t*(-368640.*t-30720.)+114944.)-57792.)+16632.)-2079.)
       /(65536.*sqrt(pow(t,7.0)*pow(t1,17.0)));
   g7=  3.0*(t*(t*(t*(t*(t*(860160.*t
         +196608.)+308480.)-177280.)+73432.)-17724.)+1899.)
       /(131072.*pow(t,4.0)*pow(t1,10.0));

  *p1= *p0*y1*(-2.0*y*g0+g1-y1*(g2-y1*(g3-y1*(g4-y1*(g5-y1*(g6-g7*y1))))));
}

/***************************************/
/*     Special Case of Airty Integral  */
/***************************************/
void omGcase2(double y, double *p0, double *p1)
{
  *p0=  1.223404016 *pow(y, 1.0/6.0)
      *(1+0.04959570165  *pow(y, -4.0/3.0)-0.008888888889 *pow(y,-2.0)
         +0.002455199181 *pow(y,-10.0/3.0)-0.0009108958061*pow(y,-4.0)
         +0.0002534684115*pow(y,-16.0/3.0));
  *p1= -0.7078817734*pow(y,-1.0/6.0)
      *(1-0.1728260369   *pow(y, -2.0/3.0)+0.0003174603174*pow(y,-2.0)
         -0.003581214850 *pow(y,- 8.0/3.0)+0.0003117824680*pow(y,-4.0)
         -0.0009073966427*pow(y,-14.0/3.0));
}

/***************************************/
/*     Asymptotic Formula              */
/***************************************/
void omGcase3(double ra, double y, double sig, double *p0, double *p1)
{
   double phi,gs,gs0,gs1,gt,gt0,gt1,ps,ps0,ps1,pt,pt0,pt1,a,b,c;
   int i,iter,flag;

   phi=ra-y*log(2*ra)+sig;
   iter=0;
   do{
       flag=0;
       gs0=gs=1.0; ps0=ps=0.0;
       gt0=gt=0.0; pt0=pt=1.0-y/ra;

       for(i=0;i<MAX_ITER;i++){
           c  =(double)(2*(i+1))*ra;
           a  =(double)(2*i+1)*y/c;
           b  =(y*y-(double)(i*(i+1)))/c;
           if((a*a+b*b)>1.0){ flag=1; break; }

           gs1=a*gs0-b*gt0;
           gt1=a*gt0+b*gs0;
           ps1=a*ps0-b*pt0-gs1/ra;
           pt1=a*pt0+b*ps0-gt1/ra;
           if(fabs(1.0-gs1/gs0)<CRIT_ITER && fabs(1.0-ps1/ps0)<CRIT_ITER)
             break;

           gs+=gs1;  gt+=gt1;
           ps+=ps1;  pt+=pt1;
           gs0=gs1;
           gt0=gt1;
           ps0=ps1;
           pt0=pt1;
       }
       if( fabs(gs*pt-ps*gt-1.0)<=CRIT_WRONSK ) break;
       ra+=5.0;
       flag=1;
       iter++;
       if(iter>MAX_ITER){
           *p0=0.0;
           *p1=0.0;
           return;
       }
   }while(flag);

  *p0=gs*cos(phi)-gt*sin(phi);
  *p1=ps*cos(phi)-pt*sin(phi);
}

/***************************************/
/*     Riccati - 2                     */
/***************************************/
void omGcase4(double ra, double y, double *p0, double *p1)
{
   double x,x1,x2,y1,y2,fai,psi,a,b,m,g0,g1,g2,g3,g4;

   x = 2.0*y/ra;
   x1= 1.0/(1.0-x);
   x2= sqrt(1.0-x);
   y1= 1.0/(2.0*y);
   y2= y1*y1;

   g0= y2*x1*x1*x1;
   g1= -(x*(9.*x-12.)+8.)/48.;
   g2= -(x*(x*(x*(x*(x*(-1890.*x+7560.)-12240.)-13440.)+16128.)-9216.)+2048.)
       /92160.;
   g3= -(x*(x*(x*(x*(x*(x*(x*(x*(x*(130977.*x
       -  873180.)+2487240.)-3588480.)+13520640.)- 9225216.)+15178880.)
       -11714560.)+6389760.)-1966080.)+  262144.)/10321920.;
   g4= 2.0*y*(x2/x+0.5*log((1.0-x2)/(1.0+x2)))+PI_4;
   fai=g4+sqrt(g0)*(g1+g0*(g2+g0*g3));

   g0*=x*x;
   g1= -(8.-3.*x)/64.;
   g2= 3.0*(x*(x*(-39.*x+208.)-448.)+1024.)/8192.;
   g3= -(x*(x*(x*(x*(-5697.*x+45576.)-159552.)+314624.)-55296.)+1105920.)
       /393216.;
   psi=x*g0*(g1+g0*(g2+g0*g3));

   g1=  (8.-3.*x)/32.;
   g2= -(x*(x*(-63.*x+336.)-704.)+1536.)/2048.;
   g3=  (x*(x*(x*(x*(-2079.*x+16632.)-57792.)+114944.)-30720.)+368640.)
       /65536.;
   a= x2/x*(1.0/x+g0*(g1+g0*(g2+g0*g3)));

   g1= -3.0*(x*(x-4.)+8.)/64.;
   g2=  3.0*(x*(x*(x*(39.*x-260.)+728.)-832.)+2560.)/4096.;
   g3= -3.0*(x*(x*(x*(x*(x*(1899.*x-17724.)+73432.)-177280.)+308480.)
           +196608.)+860160.)/131072.;
   b= x1*y1/4.0+g0*(g1+g0*(g2+g0*g3));

   m= 1.0/sqrt(x2)*exp(psi);

   *p0= m*cos(fai);
   *p1= -x*x*(a*m*sin(fai)+b*(*p0));
}


/*************************************************/
/*     N.M.Newmark                               */
/*     Proc. American Society of                 */
/*     Civil Engineering, EM 3, 67(1959)         */
/*************************************************/
void omNewmark(int n,double h, double y, double ra, double *p0, double *p1)
{
   double h2,bt,e,x,g0,g1,f0,f1,r0,r1,s1;
   int i,itr;

   y*= 2.0;
   h = -h;
   h2= h*h;
   bt= h2*BETA;
   g0= *p0;
   f0= *p1;
   r0=r1=(y/ra-1.0)*g0;
   for(i=1;i<=n;i++){
       x=ra+(double)i*h;
       itr=0;
       do{
           f1=f0+h*(r0+r1)/2.0;
           g1=g0+h*f0+h2*r0/2.0+bt*(r1-r0);
           s1=(y/x-1.0)*g1;
           e =(s1!=0.0) ? fabs(1.0-r1/s1) : fabs(s1-r1);
           itr++;
           if(itr>MAX_ITER){
               *p0=0.0;
               *p1=0.0;
               return;
           }
           r1=s1;
       }while(e>CRIT_CNVRG);
       g0=g1;
       f0=f1;
       r0=r1;
   }
   *p0=g0;
   *p1=f0;
}
