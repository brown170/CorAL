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
#include "message.h"
#include "parametermap.h"
#include "general_imager3d.h"
#include "kernel.h"
#include "kernel_chooser.h"

using namespace TNT;
using namespace std;

//--------------------------------------------
// Read/write to parameter map
//--------------------------------------------
//---------------- Read ----------------------------
bool CGeneralImager3d::Read( const parameterMap& m ){
    if (m.find("override_kernel")!=m.end()) {
        kernel_particle1 = parameter::getS(m,"kernel_particle1",kernel_particle1);
        kernel_particle2 = parameter::getS(m,"kernel_particle2",kernel_particle2);
    }
    return true;
}

//--------------- Write -----------------------------
// not used much except for diagnostics, so write everything!
bool CGeneralImager3d::Write( parameterMap& m ){
    parameter::set(m,"kernel_particle1",kernel_particle1);
    parameter::set(m,"kernel_particle2",kernel_particle2);
    return true;
}


//------------------- set_kernel -------------------------
void CGeneralImager3d::set_kernel( void ){
    parameterMap kdata;
    parameter::set(kdata,"nqmax",1000);
    parameter::set(kdata,"nrmax",100);
    parameter::set(kdata,"kellmax",4);
    parameter::set(kdata,"delq",0.5);
    parameter::set(kdata,"delr",1.0);
    parameter::set(kdata,"IDENTICAL",kernel_particle1==kernel_particle2);
    kernelPtr = chooseKernel( kernel_particle1, kernel_particle2, kdata );
}

//------------------- set_kernel -------------------------
void CGeneralImager3d::set_kernel( const parameterMap& m ){
    string kernel_particle1=parameter::getS(m,"particle1");
    string kernel_particle2=parameter::getS(m,"particle2");
    kernelPtr = chooseKernel( kernel_particle1, kernel_particle2, m );
}

//------------------- set_kernel -------------------------
void CGeneralImager3d::set_kernel( const CKernel* _kernelPtr ){
    kernelPtr = _kernelPtr;
}
