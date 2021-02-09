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
#ifndef BASISFUNCIMAGER1D_H
#define BASISFUNCIMAGER1D_H

#include <string>
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "parametermap.h"
#include "general_imager1d.h"

using namespace TNT;


//------------------------------------------
// Helper class, holds kernel matrix
//------------------------------------------
class CKernelMatrix: public Array2D<double> {

public:
    int l;
    int nq;
    double qoffset;
    double dq;
    int ncoeffs;
    bool initialized;
    string particle1;
    string particle2;

    // Constructors
    CKernelMatrix(void):l(-1),initialized(false){}

    CKernelMatrix(int _l, int _nq, double _q0, double _dq, int _nc, bool _i, string _p1, string _p2): 
        Array2D<double>(_nq,_nc,0.0),
        l(_l), 
        nq(_nq), 
        qoffset(_q0), 
        dq(_dq), 
        ncoeffs(_nc), 
        initialized(_i), 
        particle1(_p1), 
        particle2(_p2){} 

    CKernelMatrix(const CKernelMatrix& other): 
        Array2D<double>(other), 
        l(other.l), 
        nq(other.nq), 
        qoffset(other.qoffset),
        dq(other.dq), 
        ncoeffs(other.ncoeffs),  
        initialized(other.initialized), 
        particle1(other.particle1), 
        particle2(other.particle2){}
    
    // (In)Equality testers
    bool operator==(const CKernelMatrix other){
        if ((other.l<0)||(l<0)) return false;
        return 
        ( 
            l           == other.l           && 
            nq          == other.nq          && 
            qoffset     == other.qoffset     && 
            dq          == other.dq          && 
            ncoeffs     == other.ncoeffs     && 
            initialized == other.initialized && 
            particle1   == other.particle1   && 
            particle2   == other.particle2
        );
    }

    bool operator!=(const CKernelMatrix other){return !(operator==(other));}
};

//------------------------------------------
// Main imaging code
//------------------------------------------
class CBasisFuncImager1d: public CGeneralImager1d{

public:
    
    // Constructors
    CBasisFuncImager1d( void ): 
        CGeneralImager1d(),
        constrain_origin(false), 
        constrain_rmax_zero_slope(false),
        constrain_rmax_zero(false),
        kmtx(),
        conmtx(),
        convec(0),
        num_constraints(0),
        sourcePtr(0),
        __l(0), 
        __j(0){}

    // Destructors
    ~CBasisFuncImager1d( void ){}
    
    // Read/write to parameter map
    bool Read( const parameterMap& m );
    bool Write( parameterMap& m );

    // driver routines: functions that manage the (un)imaging
    bool convertCorrelationToSource( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m, const CKernel* _kernelPtr=NULL );
    bool convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn1dHisto& corrout, const parameterMap& m, const CKernel* _kernelPtr=NULL );

    // routines that do the actual (un)inversion
    virtual double imageit(CBasisFunctionExpansion1d& souout);
    virtual void unimageit(const CBasisFunctionExpansion1d& souin, CCorrFtn1dHisto& corrout);
        
    // controls for constraints
    bool constrain_origin;
    bool constrain_rmin_zero;
    bool constrain_rmax_zero_slope;
    bool constrain_rmax_zero;
    
    // internal cache of kmtx
    CKernelMatrix kmtx;

    // internal copy of constraints
    Array2D< double > conmtx;
    Array1D< double > convec;
    int num_constraints;
        
    // references to stuff needed in kp_integrand
    const CBasisFunctionExpansion1d* sourcePtr;
    int __l, __j;
    
    // functions for setting up the (un)imaging
    bool set_no_data( CSourceFtnBase& souout );
    bool initialize_source( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m );
    void set_constraints(const CBasisFunctionExpansion1d& souout);
    void set_kmtx(const CCorrFtn1dHisto& corrin, const CBasisFunctionExpansion1d& souout, const parameterMap& m );
    static int kp_integrand( unsigned int ndim, const double* x, void* classptr, unsigned fdim, double *fval);
};

#endif
