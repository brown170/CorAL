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
#include "orthogfunc_expansion1d.h"
#include <cmath>

//! Enable redimming of instances of the class.  Tries not to hose what's there.
//! Instead preserves the first min(ndata,__ndata) elements of the original data.
//! Since the different coefficients live in different subspaces, they are unaffected
//! by adding new terms.
bool COrthogonalFuncExpansion1d::setDim(int __ndata){
    ndata = __ndata;
    Array1D<double> __data(__ndata,0.0);
    Array1D<double> __uncert(__ndata,0.0);
    Array2D<double> __covmtx(__ndata,__ndata,0.0);
    double min_ndata = std::min(__ndata,data.dim());
    for (int i=0; i<min_ndata; ++i){
        for (int j=0; j<min_ndata; ++j) __covmtx[i][j] = covmtx[i][j];
        __data[i] = data[i];
        __uncert[i] = uncert[i];
    }
    data = __data;
    uncert = __uncert;
    covmtx = __covmtx;
    return true;
};

double COrthogonalFuncExpansion1d::basisFunction(double x, int i, int jderiv) const{
    if ( ( x < getLeftSupport(i) ) || (x > getRightSupport(i) ) ) return 0.0;
    double xx = rangeRemap(x);
    if (jderiv==0) return sqrt( weightFunction( xx )/normalization( i )/xScale() )*orthogFunctionValue( xx, i, 0 );
    else if (jderiv==1) {
        double w = sqrt( weightFunction( xx ) );
        return sqrt( 1.0/normalization( i )/xScale() )/xScale() * 
            ( weightFunction( xx, 1 )*orthogFunctionValue( xx, i, 0 )/2.0/w + w*orthogFunctionValue( xx, i, 1 ) );
    }
    else MESSAGE << "Did not implement basisFunction for jderiv != 0 or 1"<<ENDM_WARN;
    return 0.0;
 
}
