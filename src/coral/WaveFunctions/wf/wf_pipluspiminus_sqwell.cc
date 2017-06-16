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
#ifndef __INCLUDE_PIPLUSPIMINUS_SQWELL_CC__
#define __INCLUDE_PIPLUSPIMINUS_SQWELL_CC__
#include "wavefunction.h"

using namespace std;

CWaveFunction_pipluspiminus_sqwell::CWaveFunction_pipluspiminus_sqwell(string parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);

  m1=MPI; 
  m2=MPI; 
	IDENTICAL=0;

  q1q2=-1;

  nchannels=5;

  ellmax=2;
  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0; // S0
  ell[1]=2; // D0
  ell[2]=1; // P1
  ell[3]=0; // S2
  ell[4]=2; // D2

  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=3;
  nwells[1]=1;
  nwells[2]=3;
  nwells[3]=1;
  nwells[4]=1;

  SquareWell_MakeArrays();

  a[0][0]=0.203102; a[0][1]=0.296759; a[0][2]=1.4936;
  a[1][0]=1.1511;
  a[2][0]=0.295393; a[2][1]=0.360841; a[2][2]=0.488933;
  a[3][0]=0.188947;
  a[4][0]=0.452495;

  V0[0][0]=-28190.6; V0[0][1]=20006.6; V0[0][2]=61.6949;
  V0[1][0]=-245.695;
  V0[2][0]=-29369.6; V0[2][1]=671.979; V0[2][2]=7.1168;
  V0[3][0]=52894.6;
  V0[4][0]=16468.6;

  SquareWell_Init();

}

CWaveFunction_pipluspiminus_sqwell::~CWaveFunction_pipluspiminus_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pipluspiminus_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  const double ROOT2=sqrt(2.0);
  double psisquared;
  double P1,P2;
  double q=GetQ(iq);
  complex<double> psi,psia, psib, psisymm, psianti;
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(psia+psib)/ROOT2;
  psianti=(psia-psib)/ROOT2;
  double x;
  complex<double> DelPhi[5];

  SquareWell_GetDelPhi(iq,r,DelPhi);

	//printf("r=%g, DelPhi[2]=(%g,%g)\n",r,real(DelPhi[2]),imag(DelPhi[2]));
	
  psisquared=0.0;  
  x=q*r/HBARC;
  P1=SpherHarmonics::legendre(1,ctheta);
  P2=SpherHarmonics::legendre(2,ctheta);

  // For I=0
  psi=psisymm;
  psi+=ROOT2*DelPhi[0]/x;
  psi-=ROOT2*5.0*P2*DelPhi[1]/x;
  psisquared+=(1.0/3.0)*real(psi*conj(psi));

  // For I=1
  psi=psianti;
  psi-=ci*ROOT2*3.0*P1*DelPhi[2]/x;
  psisquared+=0.5*real(psi*conj(psi));
  
  // For I=2
  psi=psisymm;
  psi+=ROOT2*DelPhi[3]/x;
  psi-=ROOT2*5.0*P2*DelPhi[4]/x;
  psisquared+=(1.0/6.0)*real(psi*conj(psi));

	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
}

#endif
