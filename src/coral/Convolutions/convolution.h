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
#ifndef __CONVOLUTION_H
#define __CONVOLUTION_H

#include "parametermap.h"
#include "objects1d.h"
#include "corr1d_histo.h"
#include "kernel.h"
#include "crombergintegrator.h"

class CConvoluterIntegrand1d{
public:
    double q;
    int l;
    const CKernel* kernelPtr;
    const CObject1d* sourcePtr;
    CRombergIntegrator integrator;
    CConvoluterIntegrand1d(const CKernel* kPtr, const CObject1d* sPtr): 
        q(0), kernelPtr(kPtr), sourcePtr(sPtr), integrator( (void*)this, getValue )
        {l=0;if (sourcePtr) l=sourcePtr->l;}
    static double getValue(void* me, double r);
    double integrate(double rmin,double rmax);
};

CCorrFtn1dHisto convolute(const CObject1d& sourceFtn, const parameterMap& m);

#endif
