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
#include "generic_spline1d.h"
#include "linalg.h"

// Read/write to parameter map
bool CGenericSpline1d::Read(const parameterMap& m){
    vector< double > empty_vec(0);
    spline_degree = parameter::getI(m,"spline_degree",3);
    knots = stl2tntVec(parameter::getV(m,"knots_datablock",empty_vec));
    return CBasisFunctionExpansion1d::Read(m)&&(m.find("knots_datablock")!= m.end());
}

bool CGenericSpline1d::Write(parameterMap& m){
    parameter::set(m,"knots_datablock",tnt2stlVec(knots));
    parameter::set(m,"spline_degree",spline_degree);
    return CBasisFunctionExpansion1d::Write(m);
}

// CopyState
void CGenericSpline1d::CopyState(const CGenericSpline1d& a){
    CBasisFunctionExpansion1d::CopyState(a);
    knots=a.knots;
    spline_degree=a.spline_degree;
}

// Misc utility functions
int CGenericSpline1d::getKnotToLeft(double x) const{
    for (int ileft=0;ileft<knots.dim()-2;ileft++) 
        if ((knots[ileft]<x)&&(knots[ileft+1]>=x)) return ileft; 
    cerr<<"Illegal x in getKnotToLeft:"<<x<<endl;
    return -1;
}

/**
   \brief Creates default knots
   The knot list made by this routine has the following properties:
   \li the first spline_degree+1 knots are fixed to xmin
   \li the last spline_degree+1 knots are fixed to xmax
   \li the rest of the knots are equally spaced between xmin and xmax
   \warning Routine does not check that \f$ xmin <= xmax\f$
*/
bool CGenericSpline1d::setDefaultKnots( void ){
    Array1D< double > new_knots(spline_degree+1+ndata, 0.0);
    knots = new_knots;
    for (int i=0;i<spline_degree+1;++i){knots[i]=xmin;}
    for (int i=spline_degree+1;i<ndata;++i){knots[i]=xmin+(xmax-xmin)*(i-spline_degree)/(knots.dim()-2*spline_degree-1);}
    for (int i=ndata;i<knots.dim();++i){knots[i]=xmax;}
    return true;
}

/**
   \brief Creates Schoenberg-Whitney, aka optimal, knots (emulating behavior of de Boor's SPLOPT)
   \param colloc vector of collocation points to use as seed
   The knot list made by this routine has the following properties:
   \li the first bspline_degree+1 knots are fixed to colloc[0]
   \li the last bspline_degree+1 knots are fixed to colloc[colloc.dim()-1]
   \li the rest of the knots are set according to Eq. XIII.28 in de Boor's "Practical Guide to Splines"
   \warning Routine will fail if colloc not big enough & routine doesn't check if colloc is ordered like it should be 
*/
bool CGenericSpline1d::setOptimalKnots(const Array1D<double>& colloc){
    this->setDim( colloc.dim() + spline_degree - 1 );
    int first=0;
    int last=colloc.dim()-1;
    int icc=first;
    for (int i=0;i<spline_degree+1;++i){knots[i]=colloc[first];}
    for (int i=spline_degree+1;i<ndata;++i){
        knots[i]=0.;
        for (int j=0;j<spline_degree;++j){knots[i]+=colloc[icc+j]/spline_degree;}
        icc++;
    }
    for (int i=ndata;i<knots.dim();++i){knots[i]=colloc[last];}
    xmin=colloc[first];
    xmax=colloc[last];
    return true;
}


//! Tool to redim the knots, data and covmtx and keep dimensions in sync.  Will hose array contents
bool CGenericSpline1d::setDim(int ncoeffs){
    Array1D<double> new_uncert(ncoeffs,0.),new_data(ncoeffs,0.);
    Array2D<double> new_covmtx(ncoeffs,ncoeffs,0.0);
    ndata  = ncoeffs;
    data   = new_data;
    uncert = new_uncert;
    covmtx = new_covmtx;
    Array1D<double> new_knots(ncoeffs+spline_degree+1,0.);
    knots  = new_knots;
    return checkDim(spline_degree,knots.dim(),ncoeffs);
}
