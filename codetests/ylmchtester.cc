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
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

#include "constants.h"
#include "arrays.h"
#include "random.h"
#include "misc.h"
#include "parameterMap.h"

#define CHECKMULTDIV
#define CHECKXCONV

using namespace std;

int main(){
  CCHCalc chcalc;
  CRandom random(-1234);
  int lx,ly,lz,m,L,LMAX=6,i,ntest=10;
  complex<double> AfromY;
  complex<double> ci(0.0,1.0);
  double A;
  double theta,phi,x,y,z;
  
  lx=0; ly=3; lz=0;

  for(i=0;i<ntest;i++){
    theta=acos(1.0-2.0*random.ran());
    phi=2.0*PI*random.ran();

    A=chcalc.GetAFromThetaPhi(lx,ly,lz,theta,phi);
    AfromY=0.0;

    AfromY-=ci*sqrt(PI/35)*SpherHarmonics::Ylm(3,3,theta,phi);
    //AfromY-=sqrt(2*PI/105)*SpherHarmonics::Ylm(3,2,theta,phi);
    AfromY-=ci*(1.0/5.0)*sqrt(3*PI/7)*SpherHarmonics::Ylm(3,1,theta,phi);
    //AfromY-=(2.0/5.0)*sqrt(PI/7)*SpherHarmonics::Ylm(3,0,theta,phi);
    AfromY-=ci*(1.0/5.0)*sqrt(3*PI/7)*SpherHarmonics::Ylm(3,-1,theta,phi);
    //AfromY-=sqrt(2*PI/105)*SpherHarmonics::Ylm(3,-2,theta,phi);
    AfromY-=ci*sqrt(PI/35)*SpherHarmonics::Ylm(3,-3,theta,phi);

    printf("theta=%5.3f phi=%5.3f : A=%8.5f, AfromY=(%8.5f,%8.5f), error=%8.5f\n",
	   theta,phi,
	   A,real(AfromY),imag(AfromY),
	   sqrt(pow(real(AfromY)-A,2)+pow(imag(AfromY),2)));
  }

  return 0;
}
