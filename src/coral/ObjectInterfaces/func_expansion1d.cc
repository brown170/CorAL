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
#include "func_expansion1d.h"
#include <cmath>

//! Read from parameterMap
bool CBasisFunctionExpansion1d::Read(const parameterMap& m){
    xmin=parameter::getD(m,"xmin",xmin);
    xmax=parameter::getD(m,"xmax",xmax);
    return CObject1d::Read(m) && CDataSet::Read(m);
}

// write to parameterMap
bool CBasisFunctionExpansion1d::Write(parameterMap& m){
    parameter::set(m,"xmin",xmin);
    parameter::set(m,"xmax",xmax);
    return CObject1d::Write(m) && CDataSet::Write(m);
}

//! Main interface to value at a point
double CBasisFunctionExpansion1d::getValue(double x) const{
    if ( (x<xmin)||(x>xmax) ) return 0.0;
    double sum=0.0;
    for (int i=0;i<ndata;++i){sum+=data[i]*basisFunction(x,i);}
    return sum;
}

//! Main interface to uncertainty at a point
double CBasisFunctionExpansion1d::getError(double x) const{
    if ( (x<xmin)||(x>xmax) ) return 0.0;
    double sum=0.0;
    for (int i=0;i<ndata;++i){sum+=(uncert[i]*uncert[i]*basisFunction(x,i)*basisFunction(x,i));}
    return sqrt(sum);
}

//! Main interface to covariance between two points
double CBasisFunctionExpansion1d::getCovariance(double x1, double x2) const{
    if ( (x1<xmin)||(x1>xmax) || (x2<xmin)||(x2>xmax) ) return 0.0;
    double sum=0.0;
    for (int i=0;i<ndata;++i){
        if (covmtx_is_active){
            for (int j=0;j<ndata;++j) sum+=covmtx[i][j]*basisFunction(x1,i)*basisFunction(x2,j);
        }
        else sum+=uncert[i]*uncert[i]*basisFunction(x1,i)*basisFunction(x2,i);
    }
    return sum;
}

