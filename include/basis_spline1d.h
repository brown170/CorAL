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
#ifndef __BASIS_SPLINE_H__
#define __BASIS_SPLINE_H__

#include "generic_spline1d.h"

class CBasisSpline1d: public CGenericSpline1d {
public:
    //! Plain constructor, we must initialize the knots explicitly
    CBasisSpline1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=1., int bsd=3): 
        CGenericSpline1d(_l,_m,r,N,_xmin,_xmax,bsd), _splco(N,0.) 
        {Array1D<double>tmp(N+spline_degree+1,0.);knots=tmp;min_deriv=-1;max_deriv=spline_degree;}

    //! Copy constructor, copying knots taken care of by base class
    CBasisSpline1d(const CBasisSpline1d& h): 
        CGenericSpline1d(h), _splco(h._splco){}
    
    // Read/write to parameter map
    bool Read(const parameterMap& m);

    // Copy Function
    void CopyState(const CBasisSpline1d& a)
        {CGenericSpline1d::CopyState(a);_splco=a._splco;}
    
    // Main interface to source value
    double getValue(double r) const;
    double basisFunction(double x, int i, int jderiv=0) const;

    // Misc utility functions
    virtual double getLeftSupport(int i) const{return knots[i];}
    virtual double getRightSupport(int i) const{return knots[i+spline_degree+1];}

    //! Tool to redim the knots, data and covmtx and keep dimensions in sync.  Will hose array contents
    bool setDim(int ncoeffs);

protected:
    Array1D<double> _splco;
    double get_bvalue(const Array1D<double>& bcoef, double x, int jderiv) const;
    double _weight_denom(int i) const{return _weight_denom(i,spline_degree);}
    double _weight_denom(int i, int k) const{return knots[i+k]-knots[i];}
};

#endif
