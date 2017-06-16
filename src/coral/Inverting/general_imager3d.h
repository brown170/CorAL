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
#ifndef GENERALIMAGER3D_H
#define GENERALIMAGER3D_H

#include <string>
#include "parametermap.h"
#include "kernel.h"
#include "corr3d_ylm.h"
#include "soubase.h"

using namespace std;
using namespace TNT;

//------------------------------------------
// Main imaging code
//------------------------------------------
class CGeneralImager3d{

public:
    
    // Constructors
    CGeneralImager3d( void ):
        kernel_particle1(""),
        kernel_particle2(""),
        kernelPtr(NULL){}

    // Destructor
    virtual ~CGeneralImager3d( void ){}

    // Read/write to parameter map
    virtual bool Read( const parameterMap& m );
    virtual bool Write( parameterMap& m );

    // functions that actually manage the (un)imaging
    virtual bool convertCorrelationToSource( const CCorrFtn3dSphr& corrin, CSourceFtnBase& souout, const parameterMap& m )=0;
    virtual bool convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn3dSphr& corrout, const parameterMap& m )=0;

    // controls for kernel init method
    string kernel_particle1;
    string kernel_particle2;
    const CKernel* kernelPtr;
        
    //void get_usable_data( const CCorrFtn3dSphr& corrin );
    virtual bool initialize_source( const CCorrFtn3dSphr& corrin, CSourceFtnBase& souout )=0;
    virtual bool initialize_correlation( const CSourceFtnBase& souin, CCorrFtn3dSphr& corrout )=0;
    virtual bool set_no_data( CSourceFtnBase& souout )=0;
    void set_kernel( const parameterMap& m );
    void set_kernel( const CKernel* _kernelPtr );
    void set_kernel( void );
};

#endif
