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
#include "objects1d.h"
//---------------------------------------------------------
//  definitions for CObject1d
//---------------------------------------------------------

/// read from parameter map
bool CObject1d::Read(const parameterMap& s){
    l = parameter::getI(s,"l",l);
    m = parameter::getI(s,"m",m);
    realpart = parameter::getB(s,"realpart",realpart);
    realpart = !parameter::getB(s,"imagpart",!realpart);
    lx = parameter::getI(s,"lx",lx);
    ly = parameter::getI(s,"ly",ly);
    lz = parameter::getI(s,"lz",lz);
    sphrHarm = parameter::getB(s,"spherical_harmonic",sphrHarm);
    sphrHarm = !parameter::getB(s,"cartesian_harmonic",!sphrHarm);
    return true;
}

/// write to parameter map
bool CObject1d::Write(parameterMap& s){
    if (sphrHarm) {
        parameter::set(s,"spherical_harmonic",sphrHarm);
        parameter::set(s,"l",l);
        parameter::set(s,"m",m);
        parameter::set(s,"realpart",realpart);
    } else {
        parameter::set(s,"cartesian_harmonic",!sphrHarm);
        parameter::set(s,"lx",lx);
        parameter::set(s,"ly",ly);
        parameter::set(s,"lz",lz);
    }
    return true;
}

// CopyState
void CObject1d::CopyState(const CObject1d& A){
    l=A.l;
    m=A.m;
    realpart=A.realpart;
    lx=A.lx;
    ly=A.ly;
    lz=A.lz;
    sphrHarm=A.sphrHarm;
}
