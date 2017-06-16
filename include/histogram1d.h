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
#ifndef __NEW_HISTOGRAM1D_H
#define __NEW_HISTOGRAM1D_H

#include "generic_spline1d.h"
#include "parametermap.h"
#include "tnt_array1d.h"

using namespace std;

//------------------------------------------------------
//  Simple 1d binning
//------------------------------------------------------
class CHistogram1d: public CGenericSpline1d{

public:

    bool fixed_width_bins;
    
    //! Plain constructor, we must initialize the knots explicitly
    CHistogram1d(int _l=0, int _m=0, bool r=true, int N=1, double _xmin=0., double _xmax=1.): 
        CGenericSpline1d(_l,_m,r,N,_xmin,_xmax,0), fixed_width_bins(false)
        {Array1D<double>tmp(N+1,0.);knots=tmp;min_deriv=-1;max_deriv=0;}

    //! Copy constructor, copying knots taken care of by base class
    CHistogram1d(const CHistogram1d& h): 
        CGenericSpline1d(h), fixed_width_bins(h.fixed_width_bins){}

    // i/o
    bool Read(const parameterMap& m);
    bool Write(parameterMap& m);
    
    // Copy Function
    void CopyState(const CHistogram1d& a)
        {CGenericSpline1d::CopyState(a); fixed_width_bins=a.fixed_width_bins;}
    
    // Main interface to source value
    double getValue(double r) const {return data[whatBin(r)];}
    double getError(double r) const {return uncert[whatBin(r)];}
    double basisFunction(double x, int i, int jderiv=0) const;
    
    // Misc. Operations
    void setFixedWidthBins(double dx, double xoffset);
    bool inThisBin(int i, double xx) const;
    int  whatBin(double x, bool graceful=true) const;
    double binWidth(int i) const{return knots[i+1]-knots[i];}
    double leftBinEdge(int i) const{return knots[i];}
    double rightBinEdge(int i) const{return knots[i]+binWidth(i);}
    double midBin(int i) const{return knots[i]+binWidth(i)/2.;}
    
    // Inherited functions from CGenericSpline1d
    double getLeftSupport(int i) const{return leftBinEdge(i);}
    double getRightSupport(int i) const{return rightBinEdge(i);}

    //! Default knot initialization
    bool setDefaultKnots( void );
};


#endif
