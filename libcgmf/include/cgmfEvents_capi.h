/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file cgmfEvents_capi.h

  \brief CGM and CGMF C-API

*/


#ifndef __CGMFEVENTS_CAPI_H__
#define __CGMFEVENTS_CAPI_H__

#ifdef __cplusplus
extern "C" {
#endif

void setrngdptr(double (*) (void));
int  checkdatapath(const char*);
void setdatapath (const char*);

void   cgm(int, int, double, int, double []);
void   cgmf_genfissevent(int, double, double, double);
void   cgmf_genfissyields(int, double, int);

int    cgmf_getnnu();
double cgmf_getnerg(int);
double cgmf_getndircosu(int);
double cgmf_getndircosv(int);
double cgmf_getndircosw(int);
double cgmf_getntme(int);

// center-of-mass results
double cgmf_getcmnerg(int);
double cgmf_getcmndircosu(int);
double cgmf_getcmndircosv(int);
double cgmf_getcmndircosw(int);

int    cgmf_getgnu();
double cgmf_getgerg(int);
double cgmf_getgdircosu(int);
double cgmf_getgdircosv(int);
double cgmf_getgdircosw(int);
double cgmf_getgtme(int);

// light fragment
int    cgmf_getlfmass();
int    cgmf_getlfcharge();
double cgmf_getlfke();
double cgmf_getlfkepost();
double cgmf_getlfu();
float  cgmf_getlfspin();
int    cgmf_getlfparity();
int    cgmf_getlfnnu();
int    cgmf_getlfgnu();

double cgmf_getlfdircosu();
double cgmf_getlfdircosv();
double cgmf_getlfdircosw();

double cgmf_getpostlfdircosu();
double cgmf_getpostlfdircosv();
double cgmf_getpostlfdircosw();

double cgmf_getlfpremomentum_x();
double cgmf_getlfpremomentum_y();
double cgmf_getlfpremomentum_z();

double cgmf_getlfpostmomentum_x();
double cgmf_getlfpostmomentum_y();
double cgmf_getlfpostmomentum_z();

// heavy fragment

int    cgmf_gethfmass();
int    cgmf_gethfcharge();
double cgmf_gethfke();
double cgmf_gethfkepost();
double cgmf_gethfu();
float  cgmf_gethfspin();
int    cgmf_gethfparity();
int    cgmf_gethfnnu();
int    cgmf_gethfgnu();

double cgmf_gethfdircosu();
double cgmf_gethfdircosv();
double cgmf_gethfdircosw();

double cgmf_getposthfdircosu();
double cgmf_getposthfdircosv();
double cgmf_getposthfdircosw();

double cgmf_gethfpremomentum_x();
double cgmf_gethfpremomentum_y();
double cgmf_gethfpremomentum_z();

double cgmf_gethfpostmomentum_x();
double cgmf_gethfpostmomentum_y();
double cgmf_gethfpostmomentum_z();

// pre-fission neutrons

int    cgmf_getprennu();
double cgmf_getprenerg(int);
double cgmf_getprendircosu(int);
double cgmf_getprendircosv(int);
double cgmf_getprendircosw(int);

//

void cgmf_clean();

#ifdef __cplusplus
}
#endif
#endif //__CGMFEVENTS_CAPI_H__
