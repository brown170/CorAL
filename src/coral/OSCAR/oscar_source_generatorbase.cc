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
#include "oscar_source_generatorbase.h"
#include "constants.h"
#include "tnt_array1d_utils.h"

// ------------------ Read ------------------
bool COSCARSourceGeneratorBase::Read( const parameterMap& m ){
    q_cut_oscar = parameter::getD(m,"q_cut_oscar",60.0);
    return COSCARAccumulator::Read(m);
}


// ------------------ Write ------------------
bool COSCARSourceGeneratorBase::Write( parameterMap& m){
    parameter::set(m,"q_cut_oscar",q_cut_oscar);
    return COSCARAccumulator::Write(m);
}

// ------------------ pairIsGood ------------------
bool COSCARSourceGeneratorBase::pairIsGood( const COSCARLine& p1, const COSCARLine& p2 ){
    return (qinv<q_cut_oscar);
}

