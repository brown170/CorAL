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
#ifndef __INCLUDE_SOURCECALC_GX1D__
#define __INCLUDE_SOURCECALC_GX1D__
#include "sourcecalc.h"

using namespace std;

CSourceCalc_GX1D::CSourceCalc_GX1D(){
  InitSPars();
}

void CSourceCalc_GX1D::InitSPars(){
  // DEFAULT VALUES
  parameter::set(spars,"lambda",0.6);
  parameter::set(spars,"Xfrac",0.5);
  parameter::set(spars,"R",5.0);
  parameter::set(spars,"X",10.0);
  parameter::set(spars,"a",5.0);
}

void CSourceCalc_GX1D::SetSPars(double lambdaset,double Xfracset,
				double Rset,double Xset,double aset){
  parameter::set(spars,"lambda",lambdaset);
  parameter::set(spars,"R",Rset);
  parameter::set(spars,"Xfrac",Xfracset);
  parameter::set(spars,"X",Xset);
  parameter::set(spars,"a",aset);
}

void CSourceCalc_GX1D::CalcS(int lx,int ly,int lz,CCHArray *A){
  const double PI=4.0*atan(1.0);
  int ir,L;
  int nradial;
  double lambdaG,lambdaX,R,X,a,lambda,Xfrac;
  double r,z,delr,normX,normG,value;
  lambda=parameter::getD(spars,"lambda",0.3);
  Xfrac=parameter::getD(spars,"Xfrac",0.5);
  lambdaG=lambda*(1.0-Xfrac);
  lambdaX=lambda*Xfrac;
  R=parameter::getD(spars,"R",5.0);
  X=parameter::getD(spars,"X",10.0);
  a=parameter::getD(spars,"a",5.0);
  normG=pow(4.0*R*R*PI,-1.5);
  z=2.0*R*R/(X*X);
  normX=1.0/(4.0*PI*X*X*X*(2.0*z*Bessel::K1(z)+z*z*Bessel::K0(z)));

  nradial=A->GetNRADIAL();
  delr=A->GetRADSTEP();

  L=lx+ly+lz;

  for(ir=0;ir<nradial;ir++){
    r=delr*(0.5+ir);
    value=lambdaX*normX*exp(-sqrt(((r*r)/(X*X))+z*z));
    value+=lambdaG*normG*exp(-0.25*r*r/(R*R));
    if(L>0) value*=pow(r*r/(r*r+a*a),0.5*double(L));
    A->SetElement(lx,ly,lz,ir,value);
  }

}

void CSourceCalc_GX1D::CalcS(CCHArray *A){
  int dlx=1,dly=1,dlz=1;
  int lx,ly,lz,lmax;
  lmax=A->GetLMAX();
  if(A->GetXSYM()) dlx=2;
  if(A->GetYSYM()) dly=2;
  if(A->GetZSYM()) dlz=2;
  for(lx=0;lx<2;lx+=dlx){
    for(ly=0;ly<lmax-lx;ly+=dly){
      for(lz=0;lz<lmax-lx-ly;lz+=dlz){
				CalcS(lx,ly,lz,A);
      }
    }
  }
}

#endif
