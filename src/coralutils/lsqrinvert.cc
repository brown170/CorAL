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
#include <iostream>
#include "tnt_array1d.h"
#include "tnt_array2d.h"
#include "linalg.h"
#include "lsqrinvert.h"
#include "jama_qr.h"

using namespace TNT;
using namespace JAMA;
using namespace std;

#define _EPSILON 1e-8

//---------------------------------------------------------
/// Least squares w/o data covariance
//---------------------------------------------------------
int LeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double> &K,
  Array1D<double> &m, 
  Array2D<double> &Dm)
{
  if (!compare_dim(data,err,K)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=0.0;
  for (int i=0;i<data.dim();i++){
    B[i][i]=1.0/(err[i])/(err[i]);
  }
  
  // The actual algorithm
  KT=matmult(KT,B);
  Dm=svdinvert(matmult(KT,K));
  m=Dm*(KT*data);

  return 1;
}

//---------------------------------------------------------
/// Least squares w/ data covariance
//---------------------------------------------------------
int LeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double> &K,
  Array1D<double> &m, 
  Array2D<double> &Dm)
{

  if (!compare_dim(data,covmtx,K)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(covmtx);
  
  // The actual algorithm
  KT=matmult(KT,B);
  Dm=svdinvert(matmult(KT,K));
  m=Dm*(KT*data);

  return 1;
}

//---------------------------------------------------------
/// Constrained least squares w/o data covariance
/// Uses SVD to invert equations
//---------------------------------------------------------
int SVDConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm){
  
  // Get the covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=0.0;
  for (int i=0;i<data.dim();i++){
    B[i][i]=(err[i])*(err[i]);
  }
  
  return SVDConstrainedLeastSquaresInvert(data,B,K,C,c,m,Dm);
  
}

//---------------------------------------------------------
/// Constrained least squares w/ data covariance
/// Uses SVD to invert equations
//---------------------------------------------------------
int SVDConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c,
  Array1D<double>& m, 
  Array2D<double>& Dm
){

  if (!compare_dim(data,covmtx,K,C,c)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(covmtx);
  
  // Do an inversion 1st w/o constraints, so get size of model
  KT=matmult(KT,B);
  Array2D<double> Dminv(matmult(KT,K));
  Dm=svdinvert(Dminv);
  Array1D<double> x(KT*data);
  m=Dm*x;


  // Compute the chi^2's so can set tradeoff parameter
  double chi2_data=0.0, chi2_cons=0.0;
  {
    Array1D<double> dum1(K*m-data), dum2(B*dum1);
    chi2_data=inner_prod(dum1,dum2);
  }
  {
    Array1D<double> dum1(C*m-c);
    chi2_cons=inner_prod(dum1,dum1);
  }
//  cout << "Chi^2_data = "<<chi2_data<<endl;
//  cout << "Chi^2_cons = "<<chi2_cons<<endl;
  
  double trade;
   
  if ((chi2_data>_EPSILON)&&(chi2_cons>_EPSILON))
    trade=1e2*std::max(chi2_data/chi2_cons,chi2_data);
  else trade=1e2*std::max(static_cast<double>(data.dim()),chi2_data);



  // The modified algorithm
  Dm=svdinvert(Dminv+trade*matmult(transpose(C),C));
  m=Dm*(x+trade*(transpose(C)*c));
  
  
  return 1;
}

//---------------------------------------------------------
/// Constrained least squares w/o data covariance
/// Uses QR decompostion to invert equations
//---------------------------------------------------------
int QRConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K,
  const Array2D<double> &C,
  const Array1D<double> &c,
  Array1D<double> &m, 
  Array2D<double> &Dm){
  
  // Get the covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=0.0;
  for (int i=0;i<data.dim();i++){
    B[i][i]=(err[i])*(err[i]);
  }
  
  return QRConstrainedLeastSquaresInvert(data,B,K,C,c,m,Dm);
  
}

//---------------------------------------------------------
/// Constrained least squares w/ data covariance
/// Uses QR decompostion to invert equations
//---------------------------------------------------------
int QRConstrainedLeastSquaresInvert(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c,
  Array1D<double>& m, 
  Array2D<double>& Dm
){

  if (!compare_dim(data,covmtx,K,C,c)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(covmtx);
  
  // Do an inversion 1st w/o constraints, so get size of model
  KT=matmult(KT,B);
  Array2D<double> Dminv(matmult(KT,K));

  // Paste constraints onto the rest of the problem  
  Array1D<double> y(uplow_block(KT*data,c));
  Array2D<double> A(uplow_block(Dminv,C));

  // QR decomposition
  QR<double> decomp(A);
  Array2D<double> QT(transpose(decomp.getQ()));
  Array2D<double> R(decomp.getR());
  
  // JAMA's solver (does same thing as mine...)
  m=decomp.solve(y);

  // Must work a little to compute covmtx
  Dm=svdinvert(A);

  
  return !decomp.isFullRank();
}




//---------------------------------------------------------
// Dimension checkers
//---------------------------------------------------------
bool compare_dim(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K){
  // Check dimensions of everything
  if (data.dim()!=K.dim1()) 
  {
     cerr<< "data and kernel have different dimensions!"<<endl;
     return 0;
  }
  if (data.dim()!=err.dim()) 
  {
     cerr<< "data vector and error vector have different dimensions!"<<endl;
     return 0;
  }

  return 1;
}

bool compare_dim(
  const Array1D<double>& data,
  const Array1D<double>& err,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c
){

  if (!compare_dim(data,err,K)) return 0;

  // Check dimensions of everything else
  if (c.dim()!=C.dim1()) 
  {
     cerr<< "constraint vector and constraint matrix have different dimensions!"<<endl;
     return 0;
  }

  return 1;
}

bool compare_dim(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K
){  
  // Check dimensions of everything
  if (data.dim()!=K.dim1()) 
  {
     cerr<< "data and kernel have different dimensions!"<<endl;
     return 0;
  }
  if (data.dim()!=covmtx.dim1()) 
  {
     cerr<< "data vector and data covarience matrix have different dimensions!"<<endl;
     return 0;
  }
  if (covmtx.dim2()!=covmtx.dim1()) 
  {
     cerr<< "ddata covarience matrix is not square!"<<endl;
     return 0;
  }

  return 1;
}
bool compare_dim(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K,
  const Array2D<double>& C,
  const Array1D<double>& c
){  
  if (!compare_dim(data,covmtx,K)) return 0;

  // Check dimensions of everything else
  if (c.dim()!=C.dim1()) 
  {
     cerr<< "constraint vector and constraint matrix have different dimensions!"<<endl;
     return 0;
  }

  return 1;
}


//---------------------------------------------------
// CLSqrInvert Functions
//---------------------------------------------------
bool  CLSqrInvert::solve(Array1D<double> data){
    Array1D<double> fakerr(data.dim1(),1.);
    return solve(data,fakerr);
}

//-------------------------------
bool  CLSqrInvert::solve(Array1D<double> data, Array1D<double> err){
    return solve(data,err2covmtx(err));
}

//-------------------------------
bool  CLSqrInvert::solve(Array1D<double> data, Array2D<double> cov){

  if (!check_data_dim(data,cov,K)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(cov);
  
  // The actual algorithm
  KT=matmult(KT,B);
  covm=svdinvert(matmult(KT,K));
  m=covm*(KT*data);

  return true;
}


//-------------------------------
Array2D<double> CLSqrInvert::err2covmtx(Array1D<double> err){
    return makeDiagSquared(err);
}

//-------------------------------
Array2D<double>  CLSqrInvert::corrmodel(void){
    Array2D<double> corr(covm.dim1(),covm.dim2(),0.);
    Array1D<double> err(errmodel());
    for (int i=0;i<corr.dim1();i++){
        for (int j=0;j<corr.dim2();j++){
            corr[i][j]=covm[i][j]/err[i]/err[j];
        }
    }
    return corr;
}

//-------------------------------
Array1D<double>  CLSqrInvert::errmodel(void){
    Array1D<double> err(covm.dim1());
    for(int i=0;i<err.dim1();i++) err[i]=sqrt(abs(covm[i][i]));
    return err;
}

//-------------------------------
void  CLSqrInvert::usePriorModel(Array1D<double> mprior, Array2D<double> covmprior){
    if (check_prior_dim(mprior,covmprior,K)){
        useprior=true;
        priorm=mprior;
        priorcovm=covmprior;
    } else {
        cerr << "Prior model dimensions != kernel dimensions"<<endl;
        useprior=false;
    }
}

//-------------------------------
bool  CLSqrInvert::check_data_dim(
  const Array1D<double>& data,
  const Array2D<double>& covmtx,
  const Array2D<double>& K
){  
  if (data.dim()!=K.dim1()) 
  {
     cerr<< "data and kernel have different dimensions!"<<endl;
     return false;
  }
  if (data.dim()!=covmtx.dim1()) 
  {
     cerr<< "data vector and data covarience matrix have different dimensions!"<<endl;
     return false;
  }
  if (covmtx.dim2()!=covmtx.dim1()) 
  {
     cerr<< "ddata covarience matrix is not square!"<<endl;
     return false;
  }

  return true;
}

//-------------------------------
bool  CLSqrInvert::check_prior_dim(
  const Array1D<double>& mprior,
  const Array2D<double>& covmprior,
  const Array2D<double>& K
){  
  if (mprior.dim()!=K.dim2()) 
  {
     cerr<< "prior model and kernel have different dimensions!"<<endl;
     return false;
  }
  if (mprior.dim()!=covmprior.dim1()) 
  {
     cerr<< "prior model vector and covarience matrix have different dimensions!"<<endl;
     return false;
  }
  if (covmprior.dim2()!=covmprior.dim1()) 
  {
     cerr<< "prior model covarience matrix is not square!"<<endl;
     return false;
  }

  return true;
}

//---------------------------------------------------
// CLSqrInvertConstrained Functions
//---------------------------------------------------
bool  CLSqrInvertConstrained::check_constraint_dim(
    const Array2D<double>& K, 
    const Array2D<double>& conmtx,
    const Array1D<double>& convec
){
  if (conmtx.dim2()!=K.dim2()) 
  {
     cerr<< "constraint matrix and kernel have different dimensions!"<<endl;
     return false;
  }
  if (convec.dim()!=conmtx.dim1()) 
  {
     cerr<< "constraint vector and constraint matrix have different dimensions!"<<endl;
     return false;
  }

  return true;
}

//---------------------------------------------------
// CLSqrInvertSVDBigGauss Functions
//---------------------------------------------------
bool CLSqrInvertSVDBigGauss::solve(Array1D<double> data, Array2D<double> cov){

  if (!check_data_dim(data,cov,K)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(cov);
  
  // Do an inversion 1st w/o constraints, so get size of model
  KT=matmult(KT,B);
  Array2D<double> Dminv(matmult(KT,K));
  covm=svdinvert(Dminv);
  Array1D<double> x(KT*data);
  m=covm*x;


  // Compute the chi^2's so can set tradeoff parameter
  double chi2_data=0.0, chi2_cons=0.0;
  {
    Array1D<double> dum1(K*m-data), dum2(B*dum1);
    chi2_data=inner_prod(dum1,dum2);
  }
  {
    Array1D<double> dum1(conmtx*m-convec);
    chi2_cons=inner_prod(dum1,dum1);
  }
  
  double trade;
   
  if ((chi2_data>_EPSILON)&&(chi2_cons>_EPSILON))
    trade=1e2*std::max(chi2_data/chi2_cons,chi2_data);
  else trade=1e2*std::max(static_cast<double>(data.dim()),chi2_data);



  // The modified algorithm
  covm=svdinvert(Dminv+trade*matmult(transpose(conmtx),conmtx));
  m=covm*(x+trade*(transpose(conmtx)*convec));
  
  
  return 1;
}

//---------------------------------------------------
// CLSqrInvertQRBigGauss Functions
//---------------------------------------------------
bool CLSqrInvertQRBigGauss::solve(Array1D<double> data, Array2D<double> cov){

  if (!check_data_dim(data,cov,K)) exit(-1);

  // Get transpose of kernel K
  Array2D<double> KT(transpose(K));

  // Get 1/x of covarience matrix of data  
  Array2D<double> B(data.dim(),data.dim());
  B=svdinvert(cov);
  
  // Do an inversion 1st w/o constraints, so get size of model
  KT=matmult(KT,B);
  Array2D<double> Dminv(matmult(KT,K));

  // Paste constraints onto the rest of the problem  
  Array1D<double> y(uplow_block(KT*data,convec));
  Array2D<double> A(uplow_block(Dminv,conmtx));

  // QR decomposition
  QR<double> decomp(A);
  Array2D<double> QT(transpose(decomp.getQ()));
  Array2D<double> R(decomp.getR());
  
  // JAMA's solver (does same thing as mine...)
  m=decomp.solve(y);

  // Must work a little to compute covmtx
  covm=svdinvert(A);

  
  return !decomp.isFullRank();
}

//---------------------------------------------------
// CLSqrInvertSVDLagrange Functions
//---------------------------------------------------
bool CLSqrInvertSVDLagrange::solve(Array1D<double> data, Array2D<double> cov){
    if (!check_data_dim(data,cov,K)) exit(-1);

    // Get transpose of kernel K
    Array2D<double> KT(transpose(K));

    // Get 1/x of covarience matrix of data  
    Array2D<double> B(data.dim(),data.dim());
    B=svdinvert(cov);
  
    // Get KT*B, we'll use this several times
    Array2D<double> KTB(matmult(KT,B));
    
    // build "y" vector
    Array1D<double> y(uplow_block(KTB*data,convec));

    // build "M" vector
    Array2D<double> nullmtx(conmtx.dim1(),conmtx.dim1(),0.);
    Array2D<double> M(uplow_block(
        leftright_block(matmult(KTB,K), 0.5*transpose(conmtx)),
        leftright_block(        conmtx, nullmtx)
    ));
    
    // Now solve for x 
    Array2D<double> Minv(svdinvert(M));
    Array1D<double> x(Minv*y);
    
    // Unpack x 
    Array1D<double> newm(conmtx.dim2());
    Array1D<double> newlagm(x.dim()-newm.dim());
    for (int i=0;i<newm.dim();i++){newm[i]=x[i];}
    for (int i=newm.dim();i<x.dim();i++){newlagm[i-newm.dim()]=x[i];}
    m=newm;
    lagm=newlagm;
    
    // compute covariance
    Array2D<double> d2y(y.dim(),y.dim(),0.);
    Array2D<double> invcovm(matmult(KTB,K)),newcovm(m.dim(),m.dim());
    for (int i=0;i<invcovm.dim1();i++){
        for (int j=0;j<invcovm.dim2();j++) d2y[i][j]=invcovm[i][j];
    }
    Array2D<double> d2x(matmult(matmult(Minv,d2y),transpose(Minv)));
    for (int i=0;i<m.dim();i++){
        for (int j=0;j<m.dim();j++) newcovm[i][j]=d2x[i][j];
    }
    covm=newcovm;

//    covm=svdinvert(invcovm);

    return true;
}

//---------------------------------------------------
// CLSqrInvertQRLagrange Functions
//---------------------------------------------------
bool CLSqrInvertQRLagrange::solve(Array1D<double> data, Array2D<double> cov){
    if (!check_data_dim(data,cov,K)) exit(-1);
    cerr << "CLSqrInvertQRLagrange::solve not implemented yet"<<endl;
    return true;
}

