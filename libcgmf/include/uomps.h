/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file uomps.h

  \brief User-defined optical model potential library

*/

#ifndef __UOMPS_H__
#define __UOMPS_H__

#include "spline.h"

unsigned int     add_library(int,double, int, int, int, int, Optical *);

/*****************************************************/
/*   User Customized Optical Potential Parameters    */
/*****************************************************/

#define HAVE_CUSTOM_OMP 1

/*****************************************************/
/*     1) Change the number below                    */
/*****************************************************/
#define ADD_CUSTOM_OMP  4

unsigned int         scratch_omp(double, int, int, int, int, Optical *);
unsigned int             madland(double, int, int,           Optical *);
unsigned int       soukhovitskii(double, int, int, int,      Optical *);
unsigned int               spline(double, int, int, int, int, Optical *);

/*****************************************************/
/*     2) Add the definition of your function, like -*/
/*****************************************************/
/*unsigned int     your_function(double, int, int, int, int, Optical *);*/


/*****************************************************/
/*     3) Add an entry in this list                  */
/*****************************************************/
static string add_name[ADD_CUSTOM_OMP]=
  {"scratch","Madland","Efrem","spline" /*, "your_func" */};


unsigned int add_library(int id,double e,int zt,int at,int zi,int ai,Optical *omp)
{
  unsigned int potfm;
  switch(id - N_OMPLIB){
  case 1:  potfm=        scratch_omp(e,at,zt,ai,zi,omp); break;
  case 2:  potfm=            madland(e,at,zt,      omp); break;
  case 3:  potfm=      soukhovitskii(e,at,zt,   zi,omp); break;
  case 4:  potfm=             spline(e,at,zt,ai,zi,omp); break;
/*****************************************************/
/*     4) Add one line                               */
/*****************************************************/
/*case 5:  potfm=          your_func(e,at,zt,ai,zi omp); break; */
  default: potfm=0; break;

  }
  return(potfm);
}

/*****************************************************/
/*     5) And make your function                     */
/*        return value defines a potential form,     */
/*        see builtin.c                              */
/*****************************************************/
/*unsigned int your_func(double e,int at,int zt,int ai,int zi,Optical *omp)
{
  omp->v1  =  50.0 ;   omp->r0  =  1.200;   omp->rs  =  1.200;
  omp->ws1 =  10.0 ;   omp->a0  =  0.600;   omp->as  =  0.600;
   return(0x0012);
}*/




/*****************************************************/
/*         User defined OMPs                         */
/*****************************************************/
unsigned int scratch_omp(double e,int at,int zt,int ai,int zi,Optical *omp)
{
  if(ai!=1) return(0);
  if(e<0.0 || e>100.0) return(0);

  double eps = 1.0-2.0*(double)zt/(double)at;

  omp->v1  = 60.0 - 24.0*eps;
  omp->ws1 = 10.0;
  omp->vso1=  7.0;

  omp->r0  =  1.2;
  omp->rs  =  
  omp->rv =   1.2;
  omp->rvso=  1.0;

  omp->a0  =  0.6;
  omp->as  =
  omp->av  =  0.6;
  omp->avso=  0.6;

  omp->rc  =  (zi>0) ? 1.2 : 0.0;

  return(0x0013);
}


/*************************************************/
/*    Madland - Young                            */
/*        for Actinides                          */
/*************************************************/
unsigned int madland(double e, int at, int zt, Optical *omp)
{
  double eps;

  eps= 1.0-2.0*(double)zt/(double)at;

  omp->v1  = 50.378 - 27.073 *eps;
  omp->v2  = -0.354;
  omp->ws1 =  9.265 - 12.666 *eps;
  omp->ws2 = -0.232;
  omp->ws3 =  0.03318;
  omp->vso1=  6.2;

  omp->r0  =  1.264;
  omp->rs  =  1.256;
  omp->rvso=  1.100;

  omp->a0  =  0.612;
  omp->as  =  0.553 + 0.0144*e;
  omp->avso=  0.750;

  return(0x0012);
}


/*************************************************/
/*     E. Soukhovitskii S. Chiba, J.Y. Lee,      */
/*     O. Iwamoto, T. Fukahori                   */
/*     J. Nucl.Phys. G 30, 905 (2004)            */
/*************************************************/
unsigned int soukhovitskii(double e,int at, int zt, int zi, Optical *omp)
{
  double eps,ef,a13,ex,v0,v1,v2,va,vd,wsd,wsa,wvd,vso,wsod,wis,wiv,wi0,wiso,
         lmdv,lmds,lmdso,cvso,cwso,phi;
  double sn=0.0,bn=0.0;

  a13=pow((double)at,1.0/3.0);
  eps= 1.0-2.0*(double)zt/(double)at;

  if(zt == 90){
    switch(at){
    case 236:  sn = 5.11804 ; bn = 6.4381  ; break;
    case 237:  sn = 6.4381  ; bn = 4.78635 ; break;
    case 238:  sn = 4.78635 ; bn = 0.0     ; break;
    default :  sn = 0.0     ; bn = 0.0     ; break;
    }
  }else if(zt==92){
    switch(at){
    case 232:  sn = 7.250   ; bn = 5.760   ; break;
    case 233:  sn = 5.760   ; bn = 6.8437  ; break;
    case 234:  sn = 6.8437  ; bn = 5.29784 ; break;
    case 235:  sn = 5.29784 ; bn = 6.5448  ; break;
    case 236:  sn = 6.5448  ; bn = 5.1259  ; break;
    case 237:  sn = 5.1259  ; bn = 6.1520  ; break;
    case 238:  sn = 6.1520  ; bn = 4.80626 ; break;
    case 239:  sn = 4.80626 ; bn = 0.0     ; break;
    default :  sn = 0.0     ; bn = 0.0     ; break;
    }
   }else if(zt==93){
    switch(at){
    case 236:  sn = 5.740   ; bn = 6.570   ; break;
    case 237:  sn = 6.570   ; bn = 5.488   ; break;
    case 238:  sn = 5.488   ; bn = 0.0     ; break;
    default :  sn = 0.0     ; bn = 0.0     ; break;
    }
  }else if(zt==94){
    switch(at){
    case 238:  sn = 7.0005  ; bn = 5.6465  ; break;
    case 239:  sn = 5.6465  ; bn = 6.5335  ; break;
    case 240:  sn = 6.5335  ; bn = 5.2416  ; break;
    case 241:  sn = 5.2416  ; bn = 6.3094  ; break;
    case 242:  sn = 6.3094  ; bn = 0.0     ; break;
    default :  sn = 0.0     ; bn = 0.0     ; break;
    }
  }else if(zt==95){
    switch(at){
    case 241:  sn = 6.641   ; bn = 5.53757 ; break;
    case 242:  sn = 5.53757 ; bn = 0.0     ; break;
    default :  sn = 0.0     ; bn = 0.0     ; break;
    }
  }

  ef = (sn==0.0 || bn==0.0) ? -11.2814 + 0.02646 *at : -(sn+bn)*0.5;
  ex = e-ef;

  v0   =  -41.45;
  v1   =    0.03;
  v2   =    0.000205;
  vd   =   92.44;
  va   =   -0.06667;
  wsd  =   17.38;
  wsa  =    0.03833;
  wvd  =   14.74;
  vso  =    5.86;
  wsod =   -3.1;
  wis  =   11.79;
  wiv  =   81.63;
  wi0  =  100.0;
  wiso =  160.00;

  lmdv  =   0.0039075;
  lmds  =   0.01759;
  lmdso =   0.005;
  cvso  =  10.5;
  cwso  =  24.0;

/* common geometry part */
  omp->r0   = 1.245 * (1 - 0.05*(ex*ex/(ex*ex + wi0*wi0)) );
  omp->rs   = 1.2080;
  omp->rv   = 1.2476;
  omp->rvso = 1.1213;
  omp->rwso = omp->rvso;

  omp->rc   = (zi==0) ? 0.0: 1.2643;

  omp->a0   = 0.660 + 0.000253*e;
  omp->as   = 0.614;
  omp->av   = 0.594;
  omp->avso = 0.59;
  omp->awso = omp->avso;

  phi = (zi==0) ? 0.0 :
                (- v1 - 2*v2*ex + vd*lmdv*exp(-lmdv*ex))
                *(1.0 + ((zi==1)? 1: -1)*cvso*eps/(v0 + va*(at-232) + vd))
                * 0.90*(double)zt/a13;

  omp->v1   = (v0 + va*(at-232)+ v1*ex + v2*ex*ex + vd*exp(-lmdv*ex))
              *(1.0 + ((zi==1)? 1: -1)*cvso*eps/(v0 + va*(at-232) + vd)) + phi;
  omp->v2   =  0.0;
  omp->v3   =  0.0;

  omp->ws1  = (wsd + wsa*(at-232) + ((zi==1)? 1: -1)*cwso*eps)
              *exp(-lmds*ex) * ex*ex/(ex*ex+wis*wis);
  omp->ws2  = 0.0;
  omp->ws3  = 0.0;

  omp->wv1  = wvd * ex*ex/(ex*ex+wiv*wiv);
  omp->wv2  = 0.0;
  omp->wv3  = 0.0;

  omp->vso1 = vso*exp(-lmdso*ex);
  omp->vso2 = 0.0;

  omp->wso1 = wsod *ex*ex/(ex*ex+wiso*wiso);
  omp->wso2 = 0.0;

  return(0x0013);
}


#define SPLINE_FROM_FILE

/*************************************************/
/*    modified Koning-Delaroche                  */
/*    by B-spline function                       */
/*************************************************/
unsigned int spline(double e, int at, int zt, int ai, int zi, Optical *omp)
{
  unsigned int pf=koning_delaroche(e,at,zt,ai,zi,omp);

  const int ndata = 5;
  Spline spv(ndata),spw(ndata);
  double emax = 0.0;

#ifdef SPLINE_FROM_FILE
  ifstream fs("spline.dat");
  double x,y,z;

  for(int i=0 ; i<ndata ; i++){
    fs >> x;  fs >> y; fs >> z;
    spv.Setdat(x,y);
    spw.Setdat(x,z);
    if(emax<x) emax=x;
  }
  fs.close();
#else
  spv.Setdat( 0.0, 1.000);  spw.Setdat( 0.0, 1.000);
  spv.Setdat( 5.0, 1.050);  spw.Setdat( 5.0, 1.050);
  spv.Setdat(10.0, 1.100);  spw.Setdat(10.0, 1.100);
  spv.Setdat(15.0, 1.000);  spw.Setdat(15.0, 1.000);
  spv.Setdat(20.0, 0.950);  spw.Setdat(20.0, 0.950);
  emax = 20.0;
#endif

  if(e<emax){
    splineSupport(ndata,&spv);
    splineSupport(ndata,&spw);
    omp->v1  *= splineInterpolation(ndata,e,&spv);
    omp->ws1 *= splineInterpolation(ndata,e,&spw);
  }

  return(pf);
}

#endif //__UOMPS_H__
