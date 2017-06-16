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
#ifndef IMAGER_H
#define IMAGER_H

#include <string>

#include "parametermap.h"

#include "sou1d_bsplines.h"
#include "sou1d_legendre.h"
#include "sou1d_laguerre.h"
#include "sou3d_ylm.h"

#include "corr1d_histo.h"
#include "corr3d_ylm.h"

#define IS_SOURCEFUNC(s) (dynamic_cast<CSourceFtnBase*>(&s)!=NULL)

//------------------------------------------
// Main (un)imager functions, use these!  
//------------------------------------------

// --------- for 1d ----------
CSourceFtn1dBSpline image( const CCorrFtn1dHisto& corrin, const parameterMap m );
CCorrFtn1dHisto     reconstruct( const CSourceFtn1dBSpline& souin, const parameterMap m );
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, CSourceFtn1dBSpline& souout, CCorrFtn1dHisto& corrout);

CSourceFtn1dLegendrePoly image_le( const CCorrFtn1dHisto& corrin, const parameterMap m );
CCorrFtn1dHisto      reconstruct( const CSourceFtn1dLegendrePoly& souin, const parameterMap m );
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, CSourceFtn1dLegendrePoly& souout, CCorrFtn1dHisto& corrout);

CSourceFtn1dLaguerrePoly image_la( const CCorrFtn1dHisto& corrin, const parameterMap m );
CCorrFtn1dHisto      reconstruct( const CSourceFtn1dLaguerrePoly& souin, const parameterMap m );
void image_and_reconstruct( const CCorrFtn1dHisto& corrin, const parameterMap m, CSourceFtn1dLaguerrePoly& souout, CCorrFtn1dHisto& corrout);

// --------- for 3d ----------
CSourceFtn3dSphr<CSourceFtn1dBSpline>   image( const CCorrFtn3dSphr& corrin, const parameterMap m );
CCorrFtn3dSphr      reconstruct( const CSourceFtn3dSphr<CSourceFtn1dBSpline>& souin, const parameterMap m );
void image_and_reconstruct( const CCorrFtn3dSphr& corrin, const parameterMap m, CSourceFtn3dSphr<CSourceFtn1dBSpline>& souout, CCorrFtn3dSphr& corrout);


CSourceFtn3dSphr<CSourceFtn1dLegendrePoly>   image_l( const CCorrFtn3dSphr& corrin, const parameterMap m );
CCorrFtn3dSphr      reconstruct( const CSourceFtn3dSphr<CSourceFtn1dLegendrePoly>& souin, const parameterMap m );
void image_and_reconstruct( const CCorrFtn3dSphr& corrin, const parameterMap m, CSourceFtn3dSphr<CSourceFtn1dLegendrePoly>& souout, CCorrFtn3dSphr& corrout);

#endif
