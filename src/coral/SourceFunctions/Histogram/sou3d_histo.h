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
#ifndef SOURCES_3DCART_H
#define SOURCES_3DCART_H

#include <cmath>
#include "soubase.h"
#include "histogram3d.h"

class CSourceFtn3dHisto : public CSourceFtnBase, public CHistogram3d {

public:
    // create/destroy
    CSourceFtn3dHisto(string p1="", string p2="",
        int Nx=0,int Ny=0,int Nz=0,
        double Dx=1.,double Dy=1.,double Dz=1.,
        double xO=0.5, double yO=0.5, double zO=0.5): 
        CSourceFtnBase(p1,p2), CHistogram3d(Nx,Ny,Nz,Dx,Dy,Dz,xO,yO,zO){}
    CSourceFtn3dHisto(const CSourceFtn3dHisto& m) : CSourceFtnBase(m), CHistogram3d(m){}


    // read/write to parameterMap map
    bool Read(const parameterMap& s){return CSourceFtnBase::Read(s)&&CHistogram3d::Read(s);}
    bool Write(parameterMap& s){return CSourceFtnBase::Write(s)&&CHistogram3d::Write(s);}

    // copy
    void CopyState(const CSourceFtn3dHisto& m)
        {CSourceFtnBase::CopyState(m); CHistogram3d::CopyState(m);}

};

#endif
