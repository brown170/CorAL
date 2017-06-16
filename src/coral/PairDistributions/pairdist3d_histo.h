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
#ifndef PAIRDIST_3DCART_H
#define PAIRDIST_3DCART_H

#include <cmath>
#include "pairdist_base.h"
#include "histogram3d.h"

class CPairDistribution3dHisto : public CPairDistributionBase, public CHistogram3d {

public:
    // create/destroy
    CPairDistribution3dHisto(string p1="", string p2="", bool bQ=false,
        int Nx=0,int Ny=0,int Nz=0,
        double Dx=1.,double Dy=1.,double Dz=1.,
        double xO=0.5, double yO=0.5, double zO=0.5): 
        CPairDistributionBase(p1,p2,bQ), CHistogram3d(Nx,Ny,Nz,Dx,Dy,Dz,xO,yO,zO){}
    CPairDistribution3dHisto(const CPairDistribution3dHisto& m) : CPairDistributionBase(m), CHistogram3d(m){}


    // read/write to parameterMap map
    bool Read(const parameterMap& s){return CPairDistributionBase::Read(s)&&CHistogram3d::Read(s);}
    bool Write(parameterMap& s){return CPairDistributionBase::Write(s)&&CHistogram3d::Write(s);}

    // copy
    void CopyState(const CPairDistribution3dHisto& m)
        {CPairDistributionBase::CopyState(m); CHistogram3d::CopyState(m);}

};

#endif
