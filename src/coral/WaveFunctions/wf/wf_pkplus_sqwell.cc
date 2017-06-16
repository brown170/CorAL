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
#ifndef __INCLUDE_PKPLUS_SQWELL_CC__
#define __INCLUDE_PKPLUS_SQWELL_CC__
#include "wavefunction.h"
using namespace std;

CWaveFunction_pkplus_sqwell::CWaveFunction_pkplus_sqwell(string parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);

  m1=MPROTON; 
  m2=MKAON; 
	IDENTICAL=0;

  q1q2=1;
  nchannels=3;

  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0;
  ell[1]=1;
  ell[2]=1;
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=1;
  nwells[1]=2;
  nwells[2]=1;

  SquareWell_MakeArrays();

  a[0][0]=0.32857; // S11, I=1, J=1/2
  a[1][0]=0.272751; a[1][1]=1.56301; // P11, I=1, J=1/2
  a[2][0]=1.20077; // P13, I=1, J=3/2

  V0[0][0]=100912;
  V0[1][0]=679.113; V0[1][1]=41.7614;
  V0[2][0]=-45.8576;

  SquareWell_Init();
}

CWaveFunction_pkplus_sqwell::~CWaveFunction_pkplus_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pkplus_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,theta=acos(ctheta),x;
  double q=GetQ(iq);
  complex<double> psi,psia,Xlm00,Xlm10,Xlm11;
  psia=planewave[iq]->planewave(r,ctheta);
  complex<double> DelPhi[3];

  SquareWell_GetDelPhi(iq,r,DelPhi);

  x=q*r/HBARC;
  Xlm00=0.5*sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
  Xlm10=0.5*ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,0,theta,0.0)/x;
  Xlm11=0.5*ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,1,theta,0.0)/x;
  // For no spin flip
  psi=psia;
  psi+=Xlm00*DelPhi[0];
  psi+=Xlm10*((2.0/3.0)*DelPhi[2]+(1.0/3.0)*DelPhi[1]);
  psisquared=real(psi*conj(psi));
  // For spin flip
  psi=Xlm11*sqrt(2.0/9.0)*(DelPhi[2]-DelPhi[1]);
  psisquared+=real(psi*conj(psi));

	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}

#endif
