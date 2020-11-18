/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file math.h

  \brief Math definitions and data

*/


#ifndef __MATHS_H__
#define __MATHS_H__

/****************************/
/*      complex value       */
/****************************/
class Complex{ 
 public:
    double real;
    double imag;

    Complex(){
      real = 0.0;
      imag = 0.0;
    }
    Complex(double r, double i){
      real = r;
      imag = i;
    }
};


/**************************************/
/*      etc.cpp                       */
/**************************************/
Complex rational              (double, double, double, double);
Complex compsqrt              (Complex *);
Complex comppow               (Complex *, double);
double  absolute              (Complex *);


// Gauss-Legendre coefficients for numerical integration
#define MAX_GAUSS  10

static double
gauss_x[MAX_GAUSS]={
  0.9931285991850949,  0.9639719272779138,  0.9122344282513259,  0.8391169718222188,
  0.7463319064601508,  0.6360536807265150,  0.5108670019508271,  0.3737060887154196,
  0.2277858511416451,  0.0765265211334973},

gauss_a[MAX_GAUSS]={
  0.0176140071391521,  0.0406014298003869,  0.0626720483341091,  0.0832767415767047,
  0.1019301198172404,  0.1181945319615184,  0.1316886384491766,  0.1420961093183821,
  0.1491729864726037,  0.1527533871307258};


#endif //__MATH_H__
