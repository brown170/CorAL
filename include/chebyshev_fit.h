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
#ifndef CHEBYSHEV_FIT_H
#define CHEBYSHEV_FIT_H

#include "chebyshev.h"
#include "tnt_array1d.h"
#include "tnt_array2d.h"

using namespace TNT;

class CChebyshevDataFit: public CChebyshevApprox1D {

public:
    CChebyshevDataFit(void){RunQuiet();}
    ~CChebyshevDataFit(void){}
    
    /* Error on fit */
    double Error(double x) const;
    /* Set up chi2 fit*/
    void ChiSquareFitCoeffs(
        const Array1D<double> &x, 
        const Array1D<double> &y,
        const Array1D<double> &dy
        );
    void ChiSquareFitCoeffs(
        const Array1D<double> &x, 
        const Array1D<double> &y,
        const Array2D<double> &covy
        );

private:
    Array2D<double> d2c;

};


#endif
