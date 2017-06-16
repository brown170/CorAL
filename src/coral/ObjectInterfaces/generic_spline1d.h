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
#ifndef __GENERIC_SPLINE_H__
#define __GENERIC_SPLINE_H__

#include "func_expansion1d.h"
#include "tnt_array1d.h"
#include "message.h"

using namespace TNT;

class CGenericSpline1d: public CBasisFunctionExpansion1d{
public:
    Array1D<double> knots;    
    int spline_degree;

    // Constructors
    CGenericSpline1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=1.0,int deg=0): 
        CBasisFunctionExpansion1d(_l,_m,r,N,_xmin,_xmax),knots(0),spline_degree(deg){}
    CGenericSpline1d(const CGenericSpline1d& h): 
        CBasisFunctionExpansion1d(h), knots(h.knots){}
        
    // Read/write to parameter map
    bool Read(const parameterMap& m);
    bool Write(parameterMap& m);

    // CopyState
    void CopyState(const CGenericSpline1d& a);

    // Virtual functions that need to be overwritten so this class isn't virtual
    double basisFunction(double, int, int) const{
        throw MESSAGE << "Don't use CGenericSpline1d::basisFunction directly, override it in derived class!" <<ENDM_FATAL;
        return 0.0;
    }
    virtual double getLeftSupport(int i) const{
        throw MESSAGE << "Don't use CGenericSpline1d::getLeftSupport directly, override it in derived class!" <<ENDM_FATAL;
        return 0.0;
    }
    virtual double getRightSupport(int i) const{
        throw MESSAGE << "Don't use CGenericSpline1d::getRightSupport directly, override it in derived class!" <<ENDM_FATAL;
        return 0.0;
    }

    // Misc utility functions
    int getKnotToLeft(double x) const;

    //! Default knot initialization
    bool setDefaultKnots( void );

    //! Schoenberg-Whitney knots
    bool setOptimalKnots(const Array1D<double>& colloc); 

    //! Tool to redim the knots, data and covmtx and keep dimensions in sync.  Will hose array contents
    bool setDim(int ncoeffs);

    //! Check if proposed dims are legal
    bool checkDim(int splinedeg, int nknots, int ncoeffs) const
        {return (splinedeg>=0)&&(nknots>=2*splinedeg+2)&&(ncoeffs<=nknots-splinedeg-1);} 
};

#endif
