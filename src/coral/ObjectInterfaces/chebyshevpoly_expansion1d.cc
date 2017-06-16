// <<BEGIN-copyright>>
// 
//                 The GNU General Public License (GPL) Version 2, June 1991
// 
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence 
// Livermore National Laboratory. Written by Ron Soltz (soltz1@llnl.gov), David A. Brown 
// (dbrown@bnl.gov) and Scott Pratt (pratts@pa.msu.edu).
// 
// CODE-CODE-643336 All rights reserved. 
// 
// This file is part of CorAL, Version: 1.17.
// 
// Please see the file LICENSE.TXT in the main directory of this source code distribution.
// 
// This program is free software; you can redistribute it and/or modify it under the terms of 
// the GNU General Public License (as published by the Free Software Foundation) version 2, 
// dated June 1991.
// 
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
// without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
// See the terms and conditions of the GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License along with this program; 
// if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
// MA 02111-1307 USA
// 
// <<END-copyright>>
#include "chebyshevpoly_expansion1d.h"
#include <cmath>


double CChebyshevPolynomialExpansion1d::basisFunction(double x, int i, int jderiv) const{
    if ( ( x < getLeftSupport(i) ) || (x > getRightSupport(i) ) ) return 0.0;
    double xx = rangeRemap(x);
    if (jderiv==0) return sqrt( 1.0/normalization(i)/xScale() )*orthogFunctionValue(xx, i);
    if (jderiv==1) return sqrt( 1.0/normalization(i)/xScale() )/xScale()*orthogFunctionValue(xx, i, 1);
    MESSAGE << "Did not implement CChebyshevPolynomialExpansion1d::basisFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}

double CChebyshevPolynomialExpansion1d::basisFunctionInverse(double x, int i) const{
    if ( ( x < getLeftSupport(i) ) || (x > getRightSupport(i) ) ) return 0.0;
    double xx = rangeRemap(x);
    return sqrt( 1.0/normalization(i)/xScale() )*weightFunction(xx)*orthogFunctionValue(xx, i);
    MESSAGE << "Did not implement CChebyshevPolynomialExpansion1d::basisFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}


//! Main interface to value of basis function
//! Evaluates \f$ T_n(x/xscale) \f$ 
//! where \f$ T_n \f$ is the Chebyshev polynomial of order n.
double CChebyshevPolynomialExpansion1d::orthogFunctionValue(double x, int i, int jderiv) const{
    if (jderiv==0) return chebyshev_polynomial(i,x);
    if (jderiv==1) {
        if (i==0) return 0.0;
        if ( x < 1e-10 ) x=1e-10; //avoid singularity at origin
        return ( (double)i/(1-x*x) )*( x*chebyshev_polynomial(i,x) - chebyshev_polynomial(i+1,x) );
    }
    MESSAGE << "Did not implement CChebyshevPolynomialExpansion1d::orthogFunctionValue for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}
double CChebyshevPolynomialExpansion1d::weightFunction(double x, int jderiv) const{
    if (jderiv==0) return 1.0/sqrt( 1.0-x*x );
    if (jderiv==1) return x/pow( 1.0-x*x, 1.5 );
    MESSAGE << "Did not implement CChebyshevPolynomialExpansion1d::weightFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}

// Misc. functions
double CChebyshevPolynomialExpansion1d::getLeftSupport(int i) const{return xmin;}
double CChebyshevPolynomialExpansion1d::getRightSupport(int i) const{return xmax;}


double chebyshev_polynomial( int i, double x ){
    if (i < 0) return 0.0;
    double xx=x*x;
    switch (i) {
        case 0:  return 1.0;
        case 1:  return x;
        case 2:  return 2.0*xx-1.0;
        case 3:  return (4.0*xx-3.0)*x;
        case 4:  return (xx-1.0)*8.0*xx+1.0;
        case 5:  return ((16.0*xx-20.0)*xx+5.0)*x;
        case 6:  return ((32.0*xx-48.0)*xx+18.0)*xx-1.0;
        case 7:  return (((64.0*xx-112.0)*xx+56.0)*xx-7.0)*x;
        //default: return cos( static_cast<double>(i)*acos(x) );
        default: return 2.0*x*chebyshev_polynomial(i-1,x)-chebyshev_polynomial(i-2,x);
    }
    return 0.0; // Better not get here!

}
