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
#ifndef __OSCAR_CORRELATION_GENBASE_H__
#define __OSCAR_CORRELATION_GENBASE_H__

#include "parametermap.h"
#include "oscar_accumulator.h"
#include "kernel.h"
#include "random.h"

class COSCARCorrelationGeneratorBase: public COSCARAccumulator{

    public:
        
        COSCARCorrelationGeneratorBase( void ): COSCARAccumulator(), q_cut_oscar(0.060), 
            accumulation_mode("default"),ranGen(345267),
            kernelPtr(NULL), kernel_particle1(""), kernel_particle2(""){}
        ~COSCARCorrelationGeneratorBase( void ){}

        bool Read( const parameterMap& m );
        bool Write( parameterMap& m);

        bool pairIsGood( const COSCARLine& p1, const COSCARLine& p2 );
        
        void set_kernel( const parameterMap& m );
        
        double q_cut_oscar;
        
        string accumulation_mode;
        
    protected:
        CRandom ranGen;
        CKernel* kernelPtr;
        string kernel_particle1;
        string kernel_particle2;
};

#endif
