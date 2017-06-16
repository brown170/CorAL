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
#ifndef __ORTHOG_BASIS_FUNC_EXPANSION_H__
#define __ORTHOG_BASIS_FUNC_EXPANSION_H__

#include "parametermap.h"
#include "func_expansion1d.h"

class COrthogonalFuncExpansion1d: public CBasisFunctionExpansion1d {
public:

    // Constructors
    COrthogonalFuncExpansion1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=1.0): 
        CBasisFunctionExpansion1d( _l, _m, r, N, _xmin, _xmax ){}
    COrthogonalFuncExpansion1d(const COrthogonalFuncExpansion1d& h): 
        CBasisFunctionExpansion1d( h ){}

    // Destructor
    virtual ~COrthogonalFuncExpansion1d( void ){}
    
    // Read/write to parameter map
    bool Read(const parameterMap& m){return CBasisFunctionExpansion1d::Read(m);}
    bool Write(parameterMap& m){return CBasisFunctionExpansion1d::Write(m);}

    // Access to basis function of the function expansion
    double basisFunction(double x, int i, int jderiv=0) const;
    double basisFunctionInverse(double x, int i) const {return basisFunction(x,i);}  // override if not orthonormal or if need weight function
    double basisFunctionInverse(double x, int i)
        {return const_cast<const COrthogonalFuncExpansion1d*>(this)->basisFunctionInverse(x,i);}
    
    // Widgets to handle underlying orthogonal basis
    virtual double orthogFunctionValue(double x, int i, int jderiv=0) const=0;
    virtual double weightFunction(double x, int jderiv=0) const{if (jderiv==0) return 1.0; else return 0.0;}
    virtual double normalization(int j) const{return 1.0;}
    virtual double rangeRemap(double x) const{return x/xScale()-xShift();}
    virtual double xScale(void) const{return 1.0;}
    virtual double xShift(void) const{return 0.0;}

    // CopyState
    virtual void CopyState(const COrthogonalFuncExpansion1d& a)
        {CBasisFunctionExpansion1d::CopyState(a);}

    // Misc. functions
    virtual double getLeftSupport(int i) const=0;
    virtual double getRightSupport(int i) const=0;
    
    // Enable redimming of instances of the class
    bool setDim(int __ndata);

};

#endif
