/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file builtin.cpp

  \brief Built-in global optical model potentials

*/

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>

using namespace std;

#include "optical.h"
#include "etc.h"

static int  findname(string *, char *, int);


/***********************************************************/
/* Potential form bit field code description               */
/*   bit: 76543210  :  0 Imag  Volume  (WS)                */
/*                     1       Surface (deriv. WS)         */
/*                     2       Surface (gaussian)          */
/*                     3       undef.                      */
/*                     4 Real  Volume  (WS)                */
/*                     5       Surface (deriv. WS)         */
/*                     6       Surface (gaussian)          */
/*                     7       undef.                      */
/*   example                                               */
/*        00010001 = 0x0011 : Real WS + Imag WS            */
/*        00010010 = 0x0012 : Real WS + Imag dWS           */
/*        00010011 = 0x0013 : Real WS + Imag WS +dWS       */
/***********************************************************/

static const int N_OMPLIB  = 18;
/*                               Energy  At   Zt   Ai   Zi           */
static unsigned int     wilmore_hodgson(        int,      int, int, Optical *);
static unsigned int becchetti_greenlees(        int, int, int, int, Optical *);
static unsigned int            rapaport(double, int, int, int, int, Optical *);
static unsigned int         walter_guss(double, int, int, int, int, Optical *);
static unsigned int                ch89(double, int, int, int, int, Optical *);
static unsigned int    koning_delaroche(double, int, int, int, int, Optical *);
static unsigned int               perey(        int, int, int, int, Optical *);
static unsigned int    madland_schwandt(double, int, int, int, int, Optical *);
static unsigned int               lemos(                  int, int, Optical *);
static unsigned int               nolte(        int, int, int, int, Optical *);
static unsigned int           avrigeanu(double, int, int, int, int, Optical *);
static unsigned int       avrigeanu2009(double, int, int, int, int, Optical *);
static unsigned int         talys_alpha(double, int, int, int, int, Optical *);
static unsigned int            bojowald(double, int, int, int, int, Optical *);
static unsigned int           an_haixia(        int, int, int, int, Optical *);
static unsigned int           han_yinlu(        int, int, int, int, Optical *);
static unsigned int         shell_model(                       int, Optical *);
static unsigned int            test_omp(                       int, Optical *);

static string lib_name[N_OMPLIB]=
  {"Wilmore","Becchetti","Rapaport","Walter","CH89","Koning",
   "Perey","Schwandt","Lemos","Nolte","Avrigeanu","Avrigeanu2009","TALYS-A","Bojowald",
   "An", "Han","bound","test"};

unsigned int find_omp(string str)
{
   unsigned int k;
   
   k = findname(lib_name,&str[0],N_OMPLIB);
#ifdef HAVE_CUSTOM_OMP
   if(k==0){
        if( (k  = findname(add_name,&str[0],ADD_CUSTOM_OMP))!=0 ) k +=N_OMPLIB;
   }
#endif
   if(k==0)  fprintf(stderr,"ERROR     : OMP name [%s] not found\n",&str[0]);
   return(k);
}


/***********************************************************/
/*      Select Optical Potential Parameter                 */
/***********************************************************/
unsigned int omp_library(int id, int zt, int at, int zi, int ai, 
                         double e, Optical *omp)
{
   unsigned int potfm;
   switch(id){
       case  1: potfm=    wilmore_hodgson(  at,   ai,zi,omp); break;
       case  2: potfm=becchetti_greenlees(  at,zt,ai,zi,omp); break;
       case  3: potfm=           rapaport(e,at,zt,ai,zi,omp); break;
       case  4: potfm=        walter_guss(e,at,zt,ai,zi,omp); break;
       case  5: potfm=               ch89(e,at,zt,ai,zi,omp); break;
       case  6: potfm=   koning_delaroche(e,at,zt,ai,zi,omp); break;
       case  7: potfm=              perey(  at,zt,ai,zi,omp); break;
       case  8: potfm=   madland_schwandt(e,at,zt,ai,zi,omp); break;
       case  9: potfm=              lemos(        ai,zi,omp); break;
       case 10: potfm=              nolte(  at,zt,ai,zi,omp); break;
       case 11: potfm=          avrigeanu(e,at,zt,ai,zi,omp); break;
       case 12: potfm=      avrigeanu2009(e,at,zt,ai,zi,omp); break;
       case 13: potfm=        talys_alpha(e,at,zt,ai,zi,omp); break;
       case 14: potfm=           bojowald(e,at,zt,ai,zi,omp); break;
       case 15: potfm=          an_haixia(  at,zt,ai,zi,omp); break;
       case 16: potfm=          han_yinlu(  at,zt,ai,zi,omp); break;
       case 17: potfm=        shell_model(           zi,omp); break;
       case 18: potfm=           test_omp(           zi,omp); break;
       default: potfm=0                                     ; break;
   }
#ifdef HAVE_CUSTOM_OMP
   if(id>N_OMPLIB)
                potfm=     add_library(id,e,zt,at,zi,ai,omp);
#endif

   return(potfm);
}


/***********************************************************/
/*      Compare OMP Names Provided With Library            */
/***********************************************************/
int findname(string *list, char *name, int n)
{
  int  i;

  for(i=0 ; i<n ; i++){
    if(!strcmp(name,list[i].c_str())) break;
  }
  if(i==n)  return(0);
  else      return(i+1);
}


/***********************************************************/
/*      D.Wilmore, P.E.Hodgson                             */
/*      Nucl. Phys. 55, 673 (1964)                         */
/***********************************************************/
unsigned int wilmore_hodgson(int at, int ai, int zi, Optical *omp)
{
   double a;

   if(ai!=1 || zi!=0) return(0);

   a=(double)at;
   omp->v1  = 47.01 ;  omp->v2  = -0.267;  omp->v3  = -0.0018 ;
   omp->ws1 =  9.52 ;  omp->ws2 = -0.053;  omp->ws3 =  0.0    ;
   omp->r0  =  1.322-7.6e-04*a+4.0e-06*a*a-8.0e-09*a*a*a;
   omp->rs  =  1.266-3.7e-04*a+2.0e-06*a*a-4.0e-09*a*a*a;
   omp->a0  =  0.66 ;
   omp->as  =  0.48 ;
/****************  NOTE  ******************/
/* Original             Perey*            */
/* omp->v3  = -0.00118 omp->v3  = -0.0018 */
/* omp->ws2 = -0.53    omp->ws2 = -0.053  */
/*                                        */
/* *)Atomic Data and Nuclear Data Tables  */
/*   13, 293(1974)                        */
/******************************************/

   return(0x0012);
}


/***********************************************************/
/*      F.D.Becchetti, G.W.Greenlees                       */
/*      Phys. Rev. 182, 1190 (1969)                        */
/*      Ann. Rep. J.H.Williams Lab. p.116 (1969)           */
/***********************************************************/
unsigned int becchetti_greenlees(int at, int zt, int ai, int zi, Optical *omp)
{
   double eps;

   eps=1.0-2.0*(double)zt/(double)at;

   if(ai==1){ /* neutron or proton */
     if(zi==0){
       omp->v1  = 56.3-24.0*eps;
       omp->wv1 = -1.56        ;
       omp->ws1 = 13.0-12.0*eps;
       omp->rs  =  1.26 ;
       omp->rv  =  1.26 ;
       omp->rc  =  0.0  ;
       omp->as  =  0.58 ;
       omp->av  =  0.58 ;
     }
     else{
       omp->v1  = 54.0+24.0*eps +0.4*(double)zt/pow((double)at,1.0/3.0);
       omp->wv1 = -2.7         ;
       omp->ws1 = 11.8+12.0*eps;
       omp->rs  =  1.32 ;
       omp->rv  =  1.32 ;
       omp->rc  =  1.25 ;
       omp->as  =  0.51+0.7*eps;
       omp->av  =  0.51+0.7*eps;
     }
     omp->v2  = -0.32 ;
     omp->wv2 =  0.22 ;
     omp->ws2 = -0.25 ;
     omp->vso1=  6.2  ;
     omp->r0  =  1.17 ;
     omp->rvso=  1.01 ;
     omp->a0  =  0.75 ;
     omp->avso=  0.75 ;

     return(0x0013);
   }
   else if(ai==3 && zi==1){ /* triton*/
     omp->v1  = 151.9+50.0*eps;
     omp->v2  = -0.17 ;
     omp->wv1 =  41.7+44.0*eps;
     omp->wv2 = -0.33 ;
     omp->vso1=  2.5  ;

     omp->r0  =  1.20 ;
     omp->rv  =  1.40 ;
     omp->rvso=  1.20 ;
     omp->rc  =  1.3  ;

     omp->a0  =  0.72 ;
     omp->av  =  0.88 ;
     omp->avso=  0.72 ;

     return(0x0011);
   }
   else if(ai==3 && zi==2){ /* helion*/
     omp->v1  = 165.9-  6.4*eps;
     omp->v2  = -0.17 ;
     omp->wv1 =  46.0-110.0*eps;
     omp->wv2 = -0.33 ;
     omp->vso1=  2.5  ;

     omp->r0  =  1.20 ;
     omp->rv  =  1.40 ;
     omp->rvso=  1.20 ;
     omp->rc  =  1.3  ;

     omp->a0  =  0.72 ;
     omp->av  =  0.84 ;
     omp->avso=  0.72 ;

     return(0x0011);
   }
   return(0);
}


/***********************************************************/
/*      J.Rapaport, V.Kulkarni, R.W.Finlay                 */
/*      Nucl. Phys A330, 15(1979)                          */
/***********************************************************/
unsigned int rapaport(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps;

   if(ai!=1 || zi!=0) return(0);

   eps=1.0-2.0*(double)zt/(double)at;

   omp->v1  = 54.19-22.7*eps;
   omp->v2  = -0.33+0.19*eps;
   omp->v3  =  0.0;
   if(e <=15.0){
       omp->wv1 =  0.0;           omp->wv2 =  0.0 ;
       omp->ws1 =  4.28-12.8*eps; omp->ws2 =  0.4 ;
   }
   else{
       omp->wv1 = -4.3;           omp->wv2 =  0.38;
       omp->ws1 = 14.0 -10.4*eps; omp->ws2 = -0.39;
   }

   omp->vso1=  6.2  ;

   omp->r0  =  1.198;
   omp->rs  =  1.295;
   omp->rv  =  1.295;
   omp->rvso=  1.01 ;

   omp->a0  =  0.663;
   omp->as  =  0.59 ;
   omp->av  =  0.59 ;
   omp->avso=  0.75 ;

   return(0x0013);
}


/***********************************************************/
/*      R.L.Walter, P.P.Guss                               */
/*      Proc. Int. Conf. Santa Fe 1079 (1985)              */
/*      Rad. Effect, 95, 73 (1986)                         */
/***********************************************************/
unsigned int walter_guss(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,a13;

   if(ai!=1) return(0);

   eps= 1.0-2.0*(double)zt/(double)at;
   a13=pow((double)at,1.0/3.0);

   if(e<=40.0){
     if(zi==0){  omp->v1  = 52.56 -16.5  *eps;
                 omp->v2  = -0.310+ 0.081*eps;
     }
     else{       omp->v1  = 52.56 +16.5  *eps +0.4*(double)zt/a13;
                 omp->v2  = -0.310- 0.081*eps;
     }
   }
   else{
     if(zi==0)   omp->v1  = 52.56 - 0.310*40.0*(1.0+log(e/40.0))
                          -(16.50 - 0.081*40.0*(1.0+log(e/40.0)))*eps;
     else        omp->v1  = 52.56 - 0.310*40.0*(1.0+log(e/40.0))
                          +(16.50 - 0.081*40.0*(1.0+log(e/40.0)))*eps
                          +  0.4*(double)zt/a13 * 678.0/(e*e);
                 omp->v2  =  0.0;
   }

   if(e<=39.4){  omp->wv1 = -0.963;
                 omp->wv2 =  0.153;
   }
   else{         omp->wv1 = -0.963+ 0.153*e*(1.0-0.33*log(e/39.4));
                 omp->wv2 =  0.0;
   }

   if(zi==0)     omp->ws1 = 10.85 -14.94 *eps;
   else          omp->ws1 = 10.85 +14.94 *eps-1.30;
                 omp->ws2 = -0.157;

   omp->vso1=  5.767+ 2.0*eps;
   omp->vso2= -0.015;
   omp->wso1=  0.791;
   omp->wso2= -0.018;

   omp->r0  =  1.219;
   omp->rs  =  1.282; omp->rv  = 1.38 +3.76/(double)at;
   omp->rvso=  1.103; omp->rwso= 1.364;
   omp->rc  =(zi==0) ? 0.0 : 1.219;

   omp->a0  =  0.688;
   omp->as  =  0.512; omp->av  = 0.557-0.462/sqrt((double)at);
   omp->avso=  0.560; omp->awso= 0.632;
/****************  NOTE  ******************/
/* In original paper                      */
/* omp->av  = 0.557-0.462*sqrt((double)at)*/
/******************************************/

   return(0x0013);
}


/***********************************************************/
/*      R.L.Varner, et al.                                 */
/*      Phys. Rep. 201, 57(1991)                           */
/*      CH89: Chapel Hill 1989                             */
/***********************************************************/
unsigned int ch89(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,ec,a13;

   if(ai!=1) return(0);

   eps=(1.0-2.0*(double)zt/(double)at) * ((zi==0) ? -1.0 : 1.0);
   a13=pow((double)at,1.0/3.0);

   omp->v1  = 52.9  +13.1 *eps;
   omp->v2  = -0.299;
   omp->v3  =  0.0  ;
   ec       =  0.0  ;
   if(zi==1){
     omp->rc = 1.238+0.116/a13;
     ec=1.73*zt/(omp->rc*a13);
     omp->v1+= omp->v2*(-ec);
   }

   omp->wv1 =  7.8          /(1.0+exp((35.0-e+ec)/16.0));
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;
   omp->ws1 =(10.0+18.0*eps)/(1.0+exp((e-ec-36.0)/37.0));
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1=  5.9  ;

   omp->r0  =  1.250-0.225/a13;
   omp->rs  =  1.33 -0.42 /a13;
   omp->rv  =  1.33 -0.42 /a13;
   omp->rvso=  1.34 -1.2  /a13;

   omp->a0  =  0.690;
   omp->as  =  0.69 ;
   omp->av  =  0.69 ;
   omp->avso=  0.63 ;

   return(0x0013);
}


/***********************************************************/
/*      A.Koning, and J.-P.Delaroche                       */
/*      Nucl.Phys. A (2002)                                */
/*      KD2001                                             */
/***********************************************************/
unsigned int koning_delaroche(double e,int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,ef,a13,ex,v1,v2,v3,d1,d2,w1,w2,vc;

   if(ai!=1) return(0);

   a13=pow((double)at,1.0/3.0);
   eps= 1.0-2.0*(double)zt/(double)at;

/* common geometry part */
   omp->r0   = 1.3039 - 0.4054    /a13;
   omp->rs   = 1.3424 - 0.01585   *a13;
   omp->rv   = omp->r0;
   omp->rvso = 1.1854 - 0.647     /a13;
   omp->rwso = omp->rvso;
   omp->rc   = (zi==0) ? 0.0:
                         1.198 
                      +  0.697*pow((double)at, -2.0/3.0) 
                      + 12.994*pow((double)at, -5.0/3.0);

   omp->a0   = 0.6778 - 1.487e-04 *at;
   omp->av   = omp->a0;
   omp->avso = 0.59;
   omp->awso = omp->avso;

 /* for neutron */
   if(zi==0){
     ef = -11.2814             + 0.02646 *at;
     v1 = 59.30     - 21.0*eps - 0.024   *at;
     v2 = 7.228e-03            - 1.48e-06*at;
     v3 = 1.994e-05            - 2.00e-08*at;
     w1 = 12.195               + 0.0167  *at;
     d1 = 16.0      - 16.0*eps;
     vc = 0.0;

     omp->as = 0.5446 - 1.656e-04 *at;
   }
 /* for proton */
   else{
     ef = -8.4075              + 0.01378 *at;
     v1 = 59.30     + 21.0*eps - 0.024   *at;
     v2 = 7.067e-03            + 4.36e-06*at;
     v3 = 1.747e-05            + 1.50e-08*at;
     w1 = 14.336               + 0.0189  *at;
     d1 = 14.3      + 20.0*eps;
     vc = 0.42*(double)zt/a13;

     omp->as = 0.5413 + 3.963e-04 *at;
   }

   w2 = 73.55   + 0.0795*at;
   d2 =  0.0180 + 3.802e-03/(1.0+exp((at-156.0)/8.0));

   ex = e-ef;

   omp->v1  =  v1*(1.0 - ex*(v2 - ex*(v3 - ex*7.0e-09))) + vc;
   omp->v2  =  0.0;
   omp->v3  =  0.0;

   omp->ws1 =  d1*exp(-d2*ex)*ex*ex/(ex*ex+132.25);
   omp->ws2 =  0.0;
   omp->ws3 =  0.0;

   omp->wv1 =  w1            *ex*ex/(ex*ex+w2*w2);
   omp->wv2 =  0.0;
   omp->wv3 =  0.0;

   omp->vso1= (5.922 + 0.0030*at)*exp(-0.0040*ex);
   omp->vso2=  0.0;

   omp->wso1= -3.1           *ex*ex/(ex*ex+2.56e+04);
   omp->wso2=  0.0;

   return(0x0013);
}


/***********************************************************/
/*      F.G.Perey                                          */
/*      Phys. Rev. 131, 745(1963)                          */
/***********************************************************/
unsigned int perey(int at, int zt, int ai, int zi, Optical *omp)
{
   double eps,a13;

   if(ai!=1 || zi!=1) return(0);

   eps=1.0-2.0*(double)zt/(double)at;
   a13=pow((double)at,1.0/3.0);

   omp->v1  = 53.3+27.0*eps + 0.4*(double)zt/a13;
   omp->v2  = -0.55 ;
   omp->v3  =  0.0  ;

   omp->wv1 =  0.0  ;
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;

/* omp->ws1 =  3.0*a13 ; */
   omp->ws1 = 13.5  ;
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1=  7.5  ;
   omp->vso2=  0.0  ;
   omp->vso3=  0.0  ;

   omp->r0  =  1.25 ;
   omp->rs  =  1.25 ;
   omp->rv  =  0.0  ;
   omp->rvso=  1.25 ;
   omp->rc  =  1.25 ;

   omp->a0  =  0.65 ;
   omp->as  =  0.47 ;
   omp->av  =  0.0  ;
   omp->avso=  0.47 ;

   return(0x0012);
}


/***********************************************************/
/*     Madland and Schwandt                                */
/***********************************************************/
unsigned int madland_schwandt(double e, int at, int zt, int ai, int zi,
                              Optical *omp)
{
   double eps,a13,e1,e2;

   if(ai!=1 || zi!=1) return(0);

   eps=1.0-2.0*(double)zt/(double)at;
   a13=pow((double)at,1.0/3.0);
   e1 =e-80;
   e2 =log(e);

   omp->v1  =105.5*(1-0.1625*e2) + 16.5*eps;
   omp->v2  =  0.0  ;
   omp->v3  =  0.0  ;

   omp->wv1 =  6.6 +0.0273*e1 + 3.87E-06 * e1*e1*e1 ;
   omp->wv2 =  0.0  ;
   omp->wv3 =  0.0  ;

   omp->ws1 =  0.0  ;
   omp->ws2 =  0.0  ;
   omp->ws3 =  0.0  ;

   omp->vso1= 19.0*(1-0.166 *e2) - 3.75*eps;
   omp->vso2=  0.0  ;
   omp->vso3=  0.0  ;

   omp->wso1=  7.5*(1-0.248 *e2);
   omp->wso2=  0.0  ;
   omp->wso3=  0.0  ;

   omp->r0  =  1.125 + 0.001 *e;
   omp->rs  =  0.0  ;
   omp->rv  =  1.65  - 0.0024*e;
   omp->rvso=  0.92  + 0.0305*a13;
   omp->rwso=  0.877 + 0.0360*a13;
   omp->rc  =  1.25 ; /* ? */

   omp->a0  =  0.675 + 0.00031*e;
   omp->as  =  0.0  ;
   omp->av  =  0.328 + 0.00244*e;
   omp->avso=  0.768 - 0.0012 *e;
   omp->awso=  0.620;

   return(0x0011);
}


/***********************************************************/
/*      O.F.Lemos                                          */
/*      Orsay Report, Series A, No.136 (1972)              */
/***********************************************************/
unsigned int lemos(int ai, int zi, Optical *omp)
{
   if(ai!=4 || zi!=2) return(0);

   omp->v1  =193.0;    omp->v2  = -0.15;   omp->v3  =  0.0;
   omp->wv1 = 21.0;    omp->wv2 =  0.25;   omp->wv3 =  0.0;
   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;
   omp->r0  =  1.37;   omp->rv  =  1.37;   omp->rs  =  0.0;
   omp->rc  =  1.4;
   omp->a0  =  0.56;   omp->av  =  0.56;   omp->as  =  0.0;

   return(0x0011);
}


/***********************************************************/
/*      M.Nolte, H.Machner, and J.Bojowald                 */
/*      Phys. Rev. C36, 1312 (1987)                        */
/***********************************************************/
unsigned int nolte(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai!=4 || zi!=2) return(0);

   a13=pow((double)at,1.0/3.0);

   omp->v1  =101.1 + 6.051*zt/a13;    omp->v2  = -0.248;  omp->v3  =  0.0;
   omp->wv1 = 26.82- 1.706*a13   ;    omp->wv2 =  0.006;  omp->wv3 =  0.0;

   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.245;  omp->a0  =  0.817 - 0.0085*a13;
   omp->rv  =  1.570;  omp->av  =  0.692 - 0.020 *a13;
   omp->rs  =  0.0;    omp->as  =  0.0;

   omp->rc  =  1.3;
/*
    This Coulomb radius was taken from L.W.Put, Nucl.Phys. A291, 93 (1977)
*/
   return(0x0011);
}


/***********************************************************/
/*      V. Avrigeanu, P. Hodgson, M. Avrigeanu             */
/*      Phys. Rev. C49, 2136 (1994)                        */
/***********************************************************/
unsigned int avrigeanu(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai!=4 || zi!=2) return(0);

   a13=pow((double)at,1.0/3.0);

   omp->v1  =101.1 + 6.051*zt/a13;    omp->v2  = -0.248;  omp->v3  =  0.0;
   if(e>73.0){
     omp->wv1 = 26.82- 1.706*a13   ;  omp->wv2 =  0.006;  omp->wv3 =  0.0;
   }else{
     omp->wv1 = 12.64- 1.706*a13   ;  omp->wv2 =  0.200;  omp->wv3 =  0.0;
   }

   omp->ws1 =  0.0;    omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.245;  omp->a0  =  0.817 - 0.0085*a13;
   omp->rv  =  1.570;  omp->av  =  0.692 - 0.020 *a13;
   omp->rs  =  0.0;    omp->as  =  0.0;
   omp->rc  =  1.3;
/*
    This Coulomb radius was taken from L.W.Put, Nucl.Phys. A291, 93 (1977)
*/
   return(0x0011);
}


/***********************************************************/
/*      M. Avrigeanu, A.C. Obreja, F.L. Roman,             */
/*      V. Avrigeanu, W. von Oertzen                       */
/*      At. Data Nucl. Data Tables 95, 501 (2009)          */
/*                                                         */
/*      Regional OMP, 50 < A < 125, E < 50MeV              */
/***********************************************************/
unsigned int avrigeanu2009(double e, int at, int zt, int ai, int zi, Optical *omp)
{
  double a13,e1,e2,e3,rb;

  if(ai!=4 || zi!=2) return(0);

  a13=pow((double)at,1.0/3.0);

  rb = 2.66 + 1.36*a13;
  e2 = (2.59 + 10.4/(double)at)*zt/rb;
  e1 = -3.03 - 0.762*a13 + 1.24*e2;
  e3 = 23.6 + 0.181*zt/a13;


  if(e < e3){
    omp->v1  =168.0 + 0.733*zt/a13;  omp->v2  = -2.64;
  }else{
    omp->v1  =116.5 + 0.337*zt/a13;  omp->v2  = -0.453;
  }

  omp->wv1 = 2.73- 2.88*a13;  omp->wv2 =  1.11;

  if(e < e1){
   omp->ws1 =  4.0;                         omp->ws2 =  0.0;
  }else if(e < e2){
   omp->ws1 =  22.2 + 4.57*a13 - 7.446*e2;  omp->ws2 =  6.0;
  }else{
   omp->ws1 =  22.2 + 4.57*a13;             omp->ws2 = -1.446;
  }

  omp->r0  =  (e<25.0) ? 1.18 + 0.012*e : 1.48;
  omp->a0  =  (e<e2) ? 0.671 + 0.0012*at + (0.0094 - 0.0042*a13)*e2 :
                       cfmax(0.671 + 0.0012*at + (0.0094 - 0.0042*a13)*e,0.55);

  omp->rv  =  1.34;
  omp->av  =  0.50;
  omp->rs  =  1.52;
  omp->as  =  0.729 - 0.074*a13;

  omp->rc  =  1.3;

   return(0x0013);
}


/***********************************************************/
/*      A. J. Koning                                       */
/*      TALYS folding method                               */
/***********************************************************/
unsigned int talys_alpha(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   Optical nopt,popt;

   if(ai!=4 || zi!=2) return(0);

   koning_delaroche(0.25*e, at,  zt, 1, 0, &nopt);
   koning_delaroche(0.25*e, at,  zt, 1, 1, &popt);

   omp->v1   = 2*nopt.v1  + 2*popt.v1;    omp->v2  =  0.0  ;  omp->v3  =  0.0;
   omp->wv1  = 2*nopt.wv1 + 2*popt.wv1;   omp->wv2 =  0.0  ;  omp->wv3 =  0.0;
   omp->ws1  = 2*nopt.ws1 + 2*popt.ws1;   omp->ws2 =  0.0  ;  omp->ws3 =  0.0;

   omp->wso1=  0.0;    omp->wso2=  0.0;    omp->wso3=  0.0;
   omp->vso1=  0.0;    omp->vso2=  0.0;    omp->vso3=  0.0;


   omp->r0  =  (nopt.r0 + popt.r0)/2.0;
   omp->a0  =  (nopt.a0 + popt.a0)/2.0;

   omp->rv  = omp->rs  =  omp->r0;
   omp->av  = omp->as  =  omp->a0;

   omp->rc  =  popt.rc;

   return(0x0013);
}


/***********************************************************/
/*      J. Bojowald et al. for deuteron                    */
/*      Phys. Rev. C 38, 1153 (1988), Eq.(6)               */
/***********************************************************/
unsigned int bojowald(double e, int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai!=2 || zi!=1) return(0);

   a13=pow((double)at,1.0/3.0);

   omp->v1  = 81.32 + 1.43 * zt/a13;  omp->v2  = -0.240;  omp->v3  =  0.0;
   omp->wv1 = 0.132*(e-45);
   if(omp->wv1<0.0) omp->wv1 = 0.0;

   omp->ws1 =  7.80 + 1.04 * a13 - 0.712*omp->wv1;
   if(omp->ws1<0.0) omp->ws1 = 0.0;

   omp->wv2 =  0.0;    omp->wv3 =  0.0;
   omp->ws2 =  0.0;    omp->ws3 =  0.0;
   omp->vso1=  6.0;    omp->vso2=  0.0;    omp->vso3=  0.0;

   omp->r0  =  1.180;  omp->a0  =  0.636 + 0.035 *a13;
   omp->rv  =  1.270;  omp->av  =  0.768 + 0.021 *a13;

   omp->rs  =  omp->rv;
   omp->as  =  omp->av;

   omp->rvso = omp->avso = 0.78 + 0.038*a13;

   omp->rc  =  1.3;

   return(0x0013);
}


/***********************************************************/
/*      Haixia An and Chonghai Cai for deuteron            */
/*      Phys. Rev. C73 054605  (2006)                      */
/***********************************************************/
unsigned int an_haixia(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;

   if(ai!=2 || zi!=1) return(0);

   a13=pow((double)at,1.0/3.0);


   omp->v1 = 91.85 + 0.642 * zt/a13; omp->v2 = -0.249; omp->v3 = 0.000116;

   omp->wv1 = 1.104; omp->wv2 = 0.0622; omp->wv3 = 0.0;

   omp->ws1 = 10.83; omp->ws2 = -0.0306;

   omp->vso1 = 7.114; omp->vso2 = 0.0; omp->vso3 = 0.0;

   omp->r0  =  1.152 - 0.00776/a13;
   omp->a0  =  0.719 + 0.0126 *a13;

   omp->rv  =  1.305 + 0.0997/a13;
   omp->av  =  0.855 - 0.1*a13;

   omp->rs  =  1.334 + 0.152/a13;
   omp->as  =  0.531 + 0.062*a13;

   omp->rvso = 0.972;
   omp->avso = 1.011;

   omp->rc  =  1.303;

   return(0x0013);
}

/***********************************************************/
/*  Yinlu Han, Yuyang Shi and Qingbiao Shen for deuteron   */
/*      Phys. Rev. C74 044615  (2006)                      */
/***********************************************************/
unsigned int han_yinlu(int at, int zt, int ai, int zi, Optical *omp)
{
   double a13;
   double sym;

   if(ai!=2 || zi!=1) return(0);

   a13=pow((double)at,1.0/3.0);
   sym=(double)(at-2*zt)/(double)at;


   omp->v1 = 82.18 -34.811*sym + 1.058*zt/a13; omp->v2 = -0.148; omp->v3 = -0.000886;

   omp->wv1 = -4.916 + 35.0*sym; omp->wv2 = 0.0555; omp->wv3 = 0.0000442;

   omp->ws1 = 20.968 - 43.398*sym; omp->ws2 = -0.0794;

   omp->vso1 = 7.406; omp->vso2 = 0.0; omp->vso3 = 0.0;

   omp->wso1 = -0.206; omp->wso2 = 0.0; omp->wso3 = 0.0;

   omp->r0  =  1.174;
   omp->a0  =  0.809;

   omp->rv  =  1.563;
   omp->av  =  0.700 + 0.045*a13;

   omp->rs  =  1.328;
   omp->as  =  0.465 + 0.045*a13;

   omp->rvso = 1.234;
   omp->avso = 0.813;

   omp->rc  =  1.698;

   return(0x0013);
}




/***********************************************************/
/*      Shell Model Potential for Bound Nucleon            */
/*      Ross-Mark-Lawson Potential                         */
/***********************************************************/
unsigned int shell_model(int zi, Optical *omp)
{
   omp->v1  = 42.8;
   omp->vso1= 39.5/4.0;

   omp->r0  =  1.300;
   omp->rvso=  1.300;

   omp->a0  =  0.690;
   omp->avso=  0.690;
   omp->rc  = (zi==1) ? 1.25: 0.0;

   return(0x0010);
}


/***********************************************************/
/*    Test Parameters                                      */
/***********************************************************/
unsigned int test_omp(int zi, Optical *omp)
{
   omp->v1  = 50.0;   omp->v2  = -0.3;
   omp->ws1 = 10.0;   omp->ws2 =  0.4;
   omp->vso1=  7.0;

   omp->r0  =  1.200;
   omp->rs  =  1.200;
   omp->rvso=  1.200;

   omp->a0  =  0.600;
   omp->as  =  0.600;
   omp->avso=  0.600;

   omp->rc  = (zi==1) ? 1.25: 0.0;

   return(0x0012);
}
