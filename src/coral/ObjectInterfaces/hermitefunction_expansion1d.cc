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
#include "hermitefunction_expansion1d.h"
#include "constants.h"
#include <cmath>
#include "misc.h"

//! Read from parameter map
bool CHermiteFunctionExpansion1d::Read(const parameterMap& m){
    xscale = parameter::getD(m,"xscale",xscale);
    return COrthogonalFuncExpansion1d::Read(m);
}

//! Write to parameter map
bool CHermiteFunctionExpansion1d::Write(parameterMap& m){
    parameter::set(m,"xscale",xscale);
    return COrthogonalFuncExpansion1d::Write(m);
}

//! Access to basis function of the function expansion.
//! Evaluates \f$ \frac{1}{(stuff)\sqrt{xscale}} exp(-(x/xscale)^2/2.0) H_n(x/xscale) \f$ 
//! where \f$ H_n \f$ is the Hermite polynomial of order n.
//! Note, we put the weight function (the Wronskien) in here
double CHermiteFunctionExpansion1d::orthogFunctionValue(double x, int i, int jderiv) const{
    if (jderiv==0) return hermite_polynomial(i,x);
    if (jderiv==1) return 2.0*double(i)*hermite_polynomial( x, i-1 );
    MESSAGE << "Did not implement basisFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}

double CHermiteFunctionExpansion1d::weightFunction(double x, int jderiv) const{
    double y = exp(-x*x); 
    if (jderiv==0) return y; 
    if (jderiv==1) return -2.0*x*y;
    MESSAGE << "Did not implement CHermiteFunctionExpansion1d::weightFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}

double hermite_polynomial( int i, double x ){
    if (i < 0) return 0.0;
    switch (i) {
        case 0:  return 1.0;
        case 1:  return 2.0*x;
        case 2:  return 4.0*x*x-2.0;
        case 3:  return (8.0*x*x-12.0)*x;
        case 4:  return (16.0*x*x-48.0)*x*x+12.0;
        default: return 2.0*x*hermite_polynomial(i-1,x)-2.0*(i-1)*hermite_polynomial(i-2,x);
    }
    return 0.0; // Better not get here!
}

// Misc. functions
double CHermiteFunctionExpansion1d::getLeftSupport(int i) const{return xmin;}
double CHermiteFunctionExpansion1d::getRightSupport(int i) const{return xmax;}
