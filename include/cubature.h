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
#ifndef __INTEGRATECUB_H
#define __INTEGRATECUB_H

#ifndef TEST_INTEGRATOR
#include "message.h"
#endif

/* Adaptive multidimensional integration on hypercubes (or, really,
   hyper-rectangles) using cubature rules.

   A cubature rule takes a function and a hypercube and evaluates
   the function at a small number of points, returning an estimate
   of the integral as well as an estimate of the error, and also
   a suggested dimension of the hypercube to subdivide.

   Given such a rule, the adaptive integration is simple:

   1) Evaluate the cubature rule on the hypercube(s).
      Stop if converged.

   2) Pick the hypercube with the largest estimated error,
      and divide it in two along the suggested dimension.

   3) Goto (1).

*/

typedef double ( *integrand ) ( unsigned ndim, const double *x, void * );

/* Integrate the function f from xmin[dim] to xmax[dim], with at
   most maxEval function evaluations (0 for no limit),
   until the given absolute is achieved relative error.  val returns
   the integral, and estimated_error returns the estimate for the
   absolute error in val.  The return value of the function is 0
   on success and non-zero if there was an error. */
int adapt_integrate( integrand f, void *fdata,
		    unsigned dim, const double *xmin, const double *xmax, 
		    unsigned maxEval, 
		    double reqAbsError, double reqRelError, 
		    double *val, double *estimated_error );

#ifndef TEST_INTEGRATOR
//---------------------------------------
// Integrator that integrates an n-dim function
// using cubature routine
//---------------------------------------
class CIntegrateCubature{

    public:
        // Constructors and Destructors
        CIntegrateCubature(void): abserr( 1e-14 ), relerr( 1e-6 ), value( 0. ), error( 0. ), neval( 0 ),
            _ndim( 0 ), _lowerlimits( NULL ), _upperlimits( NULL ){}
        ~CIntegrateCubature(void){ delete_limits(); }

        // Initialization
        void set_ndim( int n ){ _ndim = n; new_limits(); }
        void set_limit( int n, double lolim, double uplim ){ _lowerlimits[n] = lolim; _upperlimits[n] = uplim; }

        // Member access
        int get_ndim( void ){ return _ndim; }
        double get_upper_limit( int n ){ return _upperlimits[n]; }
        double get_lower_limit( int n ){ return _lowerlimits[n]; }
        
        // Main Routine
        int compute( integrand _func, void *_funcdata ){ 
            check_ndim();
            return adapt_integrate( _func, _funcdata, _ndim, _lowerlimits, _upperlimits, neval,  abserr, relerr, &value, &error ); 
        }
        
        // Member data that the user supplies somehow
        double abserr;
        double relerr;

        // Member data we compute ourselves 
        double value;
        double error;
        int neval;

    private:
        
        // checkers for integrator parameters
        void check_ndim(void){
            if ( ( 1 > _ndim ) || ( _ndim > 15 ) ) { MESSAGE << "Bad NDIM" << ENDM_FATAL; exit( EXIT_FAILURE ); }
        }
        
        // Memory management
        void delete_limits(void){
            if ( _lowerlimits != NULL ) delete [] _lowerlimits;
            if ( _upperlimits != NULL ) delete [] _upperlimits;
            _ndim=0;
            _lowerlimits=NULL;
            _upperlimits=NULL;
        }
        void new_limits(void){
            if ( _lowerlimits != NULL ) delete_limits();
            if ( _upperlimits != NULL ) delete_limits();
            if ( _ndim <= 0 ) { MESSAGE << "Bad NDIM" << ENDM_FATAL; exit( EXIT_FAILURE ); }
            _lowerlimits = new double[_ndim];
            _upperlimits = new double[_ndim];
        }

        // Member data that the user supplies somehow
        int _ndim;
        double * _lowerlimits;
        double * _upperlimits;
};
#endif

#endif
