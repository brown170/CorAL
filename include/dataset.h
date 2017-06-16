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
#ifndef DATASET_H
#define DATASET_H

#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "parametermap.h"

using namespace TNT;

class CDataSet {
    public:
        Array1D<double> data; 
        Array1D<double> uncert;
        Array2D<double> covmtx;
        int ndata;
        bool covmtx_is_active;

        // create
        CDataSet(int N=1):
            data(N,0.),uncert(N,0.),covmtx(N,N,0.),ndata(N),covmtx_is_active(false){}
        CDataSet(const CDataSet& m): 
            data(m.data),uncert(m.uncert),covmtx(m.covmtx),ndata(m.ndata), covmtx_is_active(m.covmtx_is_active){}
        virtual ~CDataSet(void){}

        // read/write to parameter map
        bool Read(const parameterMap& com);
        bool Write(parameterMap& com);

        // misc. ftns.
        bool checkArrayDimensions(void);
        bool redim(int N);

        // error & covmtx maintainance
        void syncCovMtx(void); // syncs covariance matrix with uncertainty vector
        void syncUncert(void); // syncs uncertainty vector with covariance matrix

        // copy
        virtual void CopyState(const CDataSet& m){
            data=m.data;
            uncert=m.uncert;
            covmtx=m.covmtx;
            *const_cast<int*>(&ndata)=m.ndata;
            *const_cast<bool*>(&covmtx_is_active)=m.covmtx_is_active;
        }

};

#endif
