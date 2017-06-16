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
#ifndef SPLINEIMAGER1D_H
#define SPLINEIMAGER1D_H

#include <string>
#include "parametermap.h"
#include "basis_spline1d.h"
#include "basisfunc_imager1d.h"

#define __USE_GSL__

#ifdef __USE_GSL__
#include <gsl/gsl_multimin.h>
#endif

using namespace std;
using namespace TNT;


//----------------- zero_struct ---------------------------
//! helper struct for get_zeros function
struct zero_struct{
    double result;
    bool iszero;
};

//------------------------------------------
// Main imaging code
//------------------------------------------
class CBasisSplineImager1d: public CBasisFuncImager1d{

public:
    
    // Constructors
    CBasisSplineImager1d( void ):
        CBasisFuncImager1d(),
        qscale(50.),
        spline_degree(3),
        knot_tolerance(0.01),
        knot_init_scheme("default"),
        optimize_knots(false),
        optimize_knots_max_iter(500),
        optimize_knots_minimizer_func_weight(1.),
        optimize_knots_coellesce_knot_tol(0.01),
        optimize_knots_simplex_min_size(0.01),
        colloc_pts(0),
        user_knots(0),
        souWork(){}

    // Destructors
    ~CBasisSplineImager1d( void ){}

    // Read/write to parameter map
    bool Read( const parameterMap& m );
    bool Write( parameterMap& m );

    // routines that do the actual (un)inversion
    bool convertCorrelationToSource( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m, const CKernel* _kernelPtr=NULL );
    bool convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn1dHisto& corrout, const parameterMap& m, const CKernel* _kernelPtr=NULL );
    double imageit(CBasisFunctionExpansion1d& souout);

    // data for creating new source(s), correlation(s)
    double qscale;
    int    spline_degree;
    double knot_tolerance;

    // controls over knots
    string knot_init_scheme;
    bool   optimize_knots;
    int    optimize_knots_max_iter;
    double optimize_knots_minimizer_func_weight;
    double optimize_knots_coellesce_knot_tol;
    double optimize_knots_simplex_min_size;

    // vector for collocation points
    vector< double > colloc_pts;

    // vector of knots set by user
    vector< double > user_knots;
                
    // Disposible working source instance
    CBasisSpline1d souWork;
    
    // override this to make sure knots get initialized correctly
    bool set_no_data( CSourceFtnBase& souout );
    bool initialize_source( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m );

    // routines to set the knots
    void set_knots(CBasisSpline1d& souout);
    void set_knots_default(CBasisSpline1d& souout);
    void set_knots_from_colloc_points(CBasisSpline1d& souout);
    void set_knots_user_defined(CBasisSpline1d& souout);

    // routines to set the knots using the sampling theorem
    void set_colloc_sampling_thm(void);
    vector< double > get_zeros(double q, int nZeros, int l, double eps=1e-2);
    zero_struct find_zero(double xmin, double xmax, double q, int l, double eps);

    // routines to image and optimize knots
    double do_optimal_knots(CBasisSpline1d& souout);
#ifdef __USE_GSL__
    static double static_imageit(const gsl_vector *variableKnots, void* classPtr);
#endif
};

#endif
