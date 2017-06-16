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
#ifndef __BASIS_FUNC_EXPANSION_H__
#define __BASIS_FUNC_EXPANSION_H__

#include "objects1d.h"
#include "dataset.h"
#include "parametermap.h"

class CBasisFunctionExpansion1d: public CObject1d, public CDataSet {
public:
    double xmin, xmax;
    int min_deriv,max_deriv;

    // Constructors
    CBasisFunctionExpansion1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=100.0): 
        CObject1d(_l,_m,r), CDataSet(N), xmin(_xmin), xmax(_xmax), min_deriv(0), max_deriv(0){}
    CBasisFunctionExpansion1d(const CBasisFunctionExpansion1d& h): 
        CObject1d(h), CDataSet(h), xmin(h.xmin), xmax(h.xmax), min_deriv(h.min_deriv), max_deriv(h.max_deriv){}

    // Destructor
    virtual ~CBasisFunctionExpansion1d( void ){}
    
    // Read/write to parameter map
    bool Read(const parameterMap& m);
    bool Write(parameterMap& m);

    // Main interface to value, uncertainty and covariance at a point
    double getValue(double x) const;
    double getError(double x) const;
    double getCovariance(double x1, double x2) const;

    // Access to basis function of the function expansion
    virtual double basisFunction(double x, int i, int jderiv=0) const=0;
    double basisFunction(double x, int i, int jderiv=0)
        {return const_cast<const CBasisFunctionExpansion1d*>(this)->basisFunction(x,i,jderiv);}

    // CopyState
    virtual void CopyState(const CBasisFunctionExpansion1d& a)
        {CDataSet::CopyState(a);CObject1d::CopyState(a);xmin=a.xmin;xmax=a.xmax;}

    // Misc. functions
    virtual double getLeftSupport(int i) const=0;
    virtual double getRightSupport(int i) const=0;

    // Enable redimming of instances of the class
    virtual bool setDim(int __ndata){return CDataSet::redim(__ndata);}

    // Proposed fitting functions, since can do least squares fit with any subclass of this one
    /*
    friend bool LeastSquareFit(
        CBasisFunctionExpansion1d& A,
        const Array1D<double>& x, 
        const Array1D<double>& y, 
        const Array1D<double>& yerr);
    friend bool LeastSquareFit(
        CBasisFunctionExpansion1d& A,
        (double)func(double) );
    friend bool LeastSquareFit(
        CBasisFunctionExpansion1d& A,
        void* classPtr,
        (double)memberfunc(double) );
    */
};

#endif
