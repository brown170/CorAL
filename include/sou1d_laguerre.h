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
#ifndef SOU1D_LAGUERRE_H
#define SOU1D_LAGUERRE_H

#include "soubase.h"
#include "parametermap.h"
#include "laguerrepoly_expansion1d.h"

class CSourceFtn1dLaguerrePoly : public CLaguerrePolynomialExpansion1d, public CSourceFtnBase {
public:
    // create/destroy
    CSourceFtn1dLaguerrePoly(string p1 = "", string p2 = "", 
        int l=0, int m=0, bool repart = true, int nc=10, double rmin=0., double rmax=50., double
        rscale=5.0 ) : 
        CLaguerrePolynomialExpansion1d(l,m,repart,nc,rmin,rmax,rscale), CSourceFtnBase(p1,p2) {}
    CSourceFtn1dLaguerrePoly(const CSourceFtn1dLaguerrePoly& A) : 
        CLaguerrePolynomialExpansion1d(A), CSourceFtnBase(A) {}

    // read/write to parameter map
    bool Read(const parameterMap& s);
    bool Write(parameterMap& s);

    // copy
    void CopyState(const CSourceFtn1dLaguerrePoly& A)
        {CLaguerrePolynomialExpansion1d::CopyState(A); CSourceFtnBase::CopyState(A);}

};

#endif
