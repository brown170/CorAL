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
#include "imager.h"
#include "parametermap.h"
#include "bspline_imager1d.h"
#include "basisfunc_imager1d.h"
#include "uncoupled_imager3d.h"

//--------------------------------------------
// Main converter functions!  
//--------------------------------------------

//!-------------------- image for bsplines in 1d ------------------------
CSourceFtn1dBSpline image( const CCorrFtn1dHisto& corrin, const parameterMap m  ){
    CSourceFtn1dBSpline souout;
    CBasisSplineImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    return souout;
}


//!------------------ reconstruct for bsplines in 1d --------------------------
CCorrFtn1dHisto reconstruct( const CSourceFtn1dBSpline& souin, const parameterMap m ){
    CCorrFtn1dHisto corrout;
    CBasisSplineImager1d imager;
    imager.convertSourceToCorrelation(souin, corrout, m);
    return corrout;
}


//!------------------- image_and_reconstruct for bsplines in 1d -------------------------
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, 
    CSourceFtn1dBSpline& souout, CCorrFtn1dHisto& corrout){
    CBasisSplineImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    imager.convertSourceToCorrelation(souout, corrout, m);
}

//!-------------------- image_le for Legendre polys in 1d ------------------------
CSourceFtn1dLegendrePoly image_le( const CCorrFtn1dHisto& corrin, const parameterMap m  ){
    CSourceFtn1dLegendrePoly souout;
    CBasisFuncImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    return souout;
}


//!------------------ reconstruct for Legendre polys in 1d --------------------------
CCorrFtn1dHisto reconstruct( const CSourceFtn1dLegendrePoly& souin, const parameterMap m ){
    CCorrFtn1dHisto corrout;
    CBasisFuncImager1d imager;
    imager.convertSourceToCorrelation(souin, corrout, m);
    return corrout;
}


//!------------------- image_and_reconstruct for Legendre polys in 1d -------------------------
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, 
    CSourceFtn1dLegendrePoly& souout, CCorrFtn1dHisto& corrout){
    CBasisFuncImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    imager.convertSourceToCorrelation(souout, corrout, m);
}

//!-------------------- image_la for Laguerre polys in 1d ------------------------
CSourceFtn1dLaguerrePoly image_la( const CCorrFtn1dHisto& corrin, const parameterMap m  ){
    CSourceFtn1dLaguerrePoly souout;
    CBasisFuncImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    return souout;
}


//!------------------ reconstruct for Laguerre polys in 1d --------------------------
CCorrFtn1dHisto reconstruct( const CSourceFtn1dLaguerrePoly& souin, const parameterMap m ){
    CCorrFtn1dHisto corrout;
    CBasisFuncImager1d imager;
    imager.convertSourceToCorrelation(souin, corrout, m);
    return corrout;
}


//!------------------- image_and_reconstruct for Laguerre polys in 1d -------------------------
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, 
    CSourceFtn1dLaguerrePoly& souout, CCorrFtn1dHisto& corrout){
    CBasisFuncImager1d imager;
    imager.convertCorrelationToSource(corrin, souout, m);
    imager.convertSourceToCorrelation(souout, corrout, m);
}




//!------------------- image for bsplines in 3d -------------------------
CSourceFtn3dSphr<CSourceFtn1dBSpline> image( const CCorrFtn3dSphr& corrin, const parameterMap m  ){
    CSourceFtn3dSphr<CSourceFtn1dBSpline> souout;
    UncoupledBasisSplineImager3d imager;
    imager.convertCorrelationToSource(corrin,souout,m);
    return souout;
}

//!------------------- reconstruct for bsplines in 3d -------------------------
CCorrFtn3dSphr reconstruct( const CSourceFtn3dSphr<CSourceFtn1dBSpline>& souin, const parameterMap m ){
    CCorrFtn3dSphr corrout;
    UncoupledBasisSplineImager3d imager;
    imager.convertSourceToCorrelation(souin,corrout,m);
    return corrout;
}

//!------------------ image_and_reconstruct for bsplines in 3d --------------------------
void image_and_reconstruct( const CCorrFtn3dSphr& corrin, const parameterMap m, 
    CSourceFtn3dSphr<CSourceFtn1dBSpline>& souout, CCorrFtn3dSphr& corrout)
{
    // clear the stl::map of terms
    corrout.clear(); 
    souout.clear();
    // now image
    UncoupledBasisSplineImager3d imager;
    imager.convertCorrelationToSource(corrin,souout,m);
    imager.convertSourceToCorrelation(souout,corrout,m);
}

//!------------------- image for Legendre polys in 3d -------------------------
CSourceFtn3dSphr<CSourceFtn1dLegendrePoly> image_l( const CCorrFtn3dSphr& corrin, const parameterMap m ){
    CSourceFtn3dSphr<CSourceFtn1dLegendrePoly> souout;
    UncoupledLegendrePolyImager3d imager;
    imager.convertCorrelationToSource(corrin,souout,m);
    return souout;
}


//!------------------- reconstruct for Legendre polys in 3d -------------------------
CCorrFtn3dSphr reconstruct( const CSourceFtn3dSphr<CSourceFtn1dLegendrePoly>& souin, const parameterMap m ){
    CCorrFtn3dSphr corrout;
    UncoupledLegendrePolyImager3d imager;
    imager.convertSourceToCorrelation(souin,corrout,m);
    return corrout;
}


//!------------------ image_and_reconstruct for Legendre polys in 3d --------------------------
void image_and_reconstruct( const CCorrFtn3dSphr& corrin, const parameterMap m, CSourceFtn3dSphr<CSourceFtn1dLegendrePoly>& souout, CCorrFtn3dSphr& corrout)
{
    // clear the stl::map of terms
    corrout.clear(); 
    souout.clear();
    // now image
    UncoupledLegendrePolyImager3d imager;
    imager.convertCorrelationToSource(corrin,souout,m);
    imager.convertSourceToCorrelation(souout,corrout,m);
}
