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
#ifndef LSQRINVERT_H
#define LSQRINVERT_H
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "linalg.h"

using namespace TNT;

int LeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double> &K,
  Array1D<double> &m, 
  Array2D<double> &Dm
);
  
int LeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double> &K,
  Array1D<double> &m, 
  Array2D<double> &Dm
);

int SVDConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double> &K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm
);

int SVDConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double> &K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm
);

int QRConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double> &K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm
);

int QRConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double> &K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm
);

bool compare_dim(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K
);

bool compare_dim(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K
);

bool compare_dim(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c
);

bool compare_dim(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c
);

//---------------------------------------------------
/// Generalized Least Square Inversion
//---------------------------------------------------
class CLSqrInvert{

public:
    CLSqrInvert(Array2D<double> kern): useprior(false), K(kern), priorm(0), 
        priorcovm(1,1,0.), m(0), covm(1,1,0.){}
    virtual ~CLSqrInvert(void){}
    virtual bool solve(Array1D<double> data);
    virtual bool solve(Array1D<double> data, Array1D<double> err);
    virtual bool solve(Array1D<double> data, Array2D<double> cov);
    Array2D<double> err2covmtx(Array1D<double> err);
    Array1D<double>& model(void){return m;}
    Array2D<double>& covmodel(void){return covm;}
    Array2D<double> corrmodel(void);
    Array1D<double> errmodel(void);
    void usePriorModel(Array1D<double> mprior, Array2D<double> covmprior);
    bool check_data_dim(const Array1D<double>& data,
        const Array2D<double>& covmtx,const Array2D<double>& K);
    bool check_prior_dim(const Array1D<double>& mprior,
        const Array2D<double>& covmprior,const Array2D<double>& K);

protected :
    bool useprior;
    Array2D<double> K;
    Array1D<double> priorm;
    Array2D<double> priorcovm;
    Array1D<double> m;
    Array2D<double> covm;

};

//---------------------------------------------------
/// Base class for ...
/// Generalized Least Square Inversion w/ Constraints
//---------------------------------------------------
class CLSqrInvertConstrained: public CLSqrInvert{

public:
    CLSqrInvertConstrained(Array2D<double> kern, Array2D<double> cmtx,
        Array1D<double> cvec): CLSqrInvert(kern), conmtx(cmtx), convec(cvec){
        if (!check_constraint_dim(kern,cmtx,cvec)) exit(-1);
    }
    virtual ~CLSqrInvertConstrained(void){}
    virtual bool solve(Array1D<double> data, Array2D<double> cov)=0;
    bool check_constraint_dim(const Array2D<double>& K, 
        const Array2D<double>& conmtx,const Array1D<double>& convec);
        
protected:
    Array2D<double> conmtx;
    Array1D<double> convec;
   
};

//---------------------------------------------------
/// Generalized Least Square Inversion w/ Constraints
/// using delta-ftn -> Gaussian trick and SVD decomp
//---------------------------------------------------
class CLSqrInvertSVDBigGauss: public CLSqrInvertConstrained{

public:
    CLSqrInvertSVDBigGauss(Array2D<double> kern, Array2D<double> cmtx,
        Array1D<double> cvec): CLSqrInvertConstrained(kern,cmtx,cvec){}
    ~CLSqrInvertSVDBigGauss(void){}
    bool solve(Array1D<double> data, Array2D<double> cov);

};

//---------------------------------------------------
/// Generalized Least Square Inversion w/ Constraints
/// using delta-ftn -> Gaussian trick and QR decomp
//---------------------------------------------------
class CLSqrInvertQRBigGauss: public CLSqrInvertConstrained{

public:
    CLSqrInvertQRBigGauss(Array2D<double> kern, Array2D<double> cmtx,
        Array1D<double> cvec): CLSqrInvertConstrained(kern,cmtx,cvec){}
    ~CLSqrInvertQRBigGauss(void){}
    bool solve(Array1D<double> data, Array2D<double> cov);

};

//---------------------------------------------------
/// Generalized Least Square Inversion w/ Constraints
/// using Lagrange multiplier trick and SVD decomp
//---------------------------------------------------
class CLSqrInvertSVDLagrange: public CLSqrInvertConstrained{

public:
    CLSqrInvertSVDLagrange(Array2D<double> kern, Array2D<double> cmtx,
        Array1D<double> cvec): CLSqrInvertConstrained(kern,cmtx,cvec),lagm(0){}
    ~CLSqrInvertSVDLagrange(void){}
    bool solve(Array1D<double> data, Array2D<double> cov);
    Array1D<double> lagrange_multipliers(void){return lagm;}
protected:
    Array1D<double> lagm;
};

//---------------------------------------------------
/// Generalized Least Square Inversion w/ Constraints
/// using Lagrange multiplier trick and QR decomp
//---------------------------------------------------
class CLSqrInvertQRLagrange: public CLSqrInvertConstrained{

public:
    CLSqrInvertQRLagrange(Array2D<double> kern, Array2D<double> cmtx,
        Array1D<double> cvec): CLSqrInvertConstrained(kern,cmtx,cvec),lagm(0){}
    ~CLSqrInvertQRLagrange(void){}
    bool solve(Array1D<double> data, Array2D<double> cov);
    Array1D<double> lagrange_multipliers(void){return lagm;}
protected:
    Array1D<double> lagm;

};

#endif
