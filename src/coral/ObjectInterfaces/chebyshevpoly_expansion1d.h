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
#ifndef __CHEBYSHEV_POLY_EXPANSION_H__
#define __CHEBYSHEV_POLY_EXPANSION_H__

#include "parametermap.h"
#include "orthogfunc_expansion1d.h"
#include "constants.h"

double chebyshev_polynomial( int i, double x );

class CChebyshevPolynomialExpansion1d: public COrthogonalFuncExpansion1d {
public:

    // Constructors
    CChebyshevPolynomialExpansion1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=1.0): 
        COrthogonalFuncExpansion1d( _l, _m, r, N, _xmin, _xmax ){}
    CChebyshevPolynomialExpansion1d(const CChebyshevPolynomialExpansion1d& h): 
        COrthogonalFuncExpansion1d( h ){}

    // Destructor
    ~CChebyshevPolynomialExpansion1d( void ){}
    
    // Read/write to parameter map
    bool Read(const parameterMap& m){return COrthogonalFuncExpansion1d::Read(m);}
    bool Write(parameterMap& m){return COrthogonalFuncExpansion1d::Write(m);}

    double basisFunction(double x, int i, int jderiv=0) const;
    double basisFunctionInverse(double x, int i) const;  // override if not orthonormal or if need weight function
    double orthogFunctionValue(double x, int i, int jderiv=0) const;
    double weightFunction(double x, int jderiv=0) const;
    double normalization(int j) const{if (j==0) return PI; else return PI/2.0;}
    double xScale(void) const{return (xmax-xmin)/2.0;}
    double xShift(void) const{return (xmax+xmin)/(xmax-xmin);}

    // CopyState
    void CopyState(const CChebyshevPolynomialExpansion1d& a)
        {COrthogonalFuncExpansion1d::CopyState(a);}

    // Misc. functions
    double getLeftSupport(int i) const;
    double getRightSupport(int i) const;

};

#endif
