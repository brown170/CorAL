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
#include "legendrepoly_expansion1d.h"
#include <gsl/gsl_sf_legendre.h>


//! Main interface to value of basis function
//! Evaluates \f$ L_n(x/xscale) \f$ 
//! where \f$ L_n \f$ is the Legendre polynomial of order n.
double CLegendrePolynomialExpansion1d::orthogFunctionValue(double x, int i, int jderiv) const{
    if (jderiv==0) return gsl_sf_legendre_Pl(i,x);
    if (jderiv==1) {
        if (i==0) return 0.0;
        if ( x < 1e-10 ) x=1e-10; //avoid singularity at origin
        return ( (double)i/(1-x*x) )*( gsl_sf_legendre_Pl(i-1,x) - x*gsl_sf_legendre_Pl(i,x) );
    }
    MESSAGE << "Did not implement basisFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
}


// Misc. functions
double CLegendrePolynomialExpansion1d::getLeftSupport(int i) const{return xmin;}
double CLegendrePolynomialExpansion1d::getRightSupport(int i) const{return xmax;}
