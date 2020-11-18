/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file gdr.h

  \brief Prototype of functions for global GDR parameters

*/
#ifndef __GDR_H__
#define __GDR_H__

#include <string>

using namespace std;

class GDR{
 private:
    string     XL            ;     /* E or M, and multipolarity        */
    double     energy        ;     /* GDR energy                       */
    double     width         ;     /* GDR width                        */
    double     sigma         ;     /* GDR peak-cross section           */
 public:
    GDR(){
      XL      = "  ";
      energy  = 0.0;
      width   = 0.0;
      sigma   = 0.0;
    }
    void setGDR(string em, double e, double w, double s){
      XL      = em;
      energy  = e;
      width   = w;
      sigma   = s;
    }
    void clear(){
      XL      = "  ";
      energy  = 0.0;
      width   = 0.0;
      sigma   = 0.0;
    }
    string getXL    () {return XL;}
    double getEnergy() {return energy;}
    double getWidth () {return width ;}
    double getSigma () {return sigma ;}
    char   getEM    () {char em = (XL[0]=='E') ? 'E' : 'M'; return(em);}
    int    getL     () {return( (int)(XL[1])-'0');}
    void   setXL    (string em) {XL     = em;}
    void   setEnergy(double e ) {energy = e;}
    void   setWidth (double w ) {width  = w;}
    void   setSigma (double s ) {sigma  = s;}
};


/**************************************/
/*      gdr.cpp                       */
/**************************************/
void    gdrDArigo              (double, double, GDR *);
void    gdrDArigoDoubleHump0   (double, double, GDR *);
void    gdrDArigoDoubleHump1   (double, double, GDR *);
void    gdrM1                  (double, GDR *);
void    gdrE2                  (double, double, GDR *);

#endif //__GDR_H__
