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
#include "oscar_correlation_generatorbase.h"
#include "kernel_chooser.h"

// ------------------ set_kernel ------------------
void COSCARCorrelationGeneratorBase::set_kernel( const parameterMap& m ){
    kernel_particle1 = parameter::getS(m,"particle1",kernel_particle1);
    kernel_particle2 = parameter::getS(m,"particle2",kernel_particle2);
    kernelPtr = chooseKernel( kernel_particle1, kernel_particle2, m );
}

// ------------------ Read ------------------
bool COSCARCorrelationGeneratorBase::Read( const parameterMap& m ){
    accumulation_mode = parameter::getS(m,"correlation_accumulation_mode",accumulation_mode);
    q_cut_oscar = parameter::getD(m,"q_cut_oscar",60.0);
    return COSCARAccumulator::Read(m);
}

// ------------------ Write ------------------
bool COSCARCorrelationGeneratorBase::Write( parameterMap& m){
    parameter::set(m,"correlation_accumulation_mode",accumulation_mode);
    parameter::set(m,"q_cut_oscar",q_cut_oscar);
    return COSCARAccumulator::Write(m);
}

// ------------------ pairIsGood ------------------
bool COSCARCorrelationGeneratorBase::pairIsGood( const COSCARLine& p1, const COSCARLine& p2 ){
    return (qinv<q_cut_oscar);
}
