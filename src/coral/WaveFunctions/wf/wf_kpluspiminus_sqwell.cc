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
#ifndef __INCLUDE_KPLUSPIMINUS_SQWELL__
#define __INCLUDE_KPLUSPIMINUS_SQWELL__
#include "wavefunction.h"

using namespace std;

CWaveFunction_kpluspiminus_sqwell::CWaveFunction_kpluspiminus_sqwell(string parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);

  m1=MPI; 
  m2=MKAON;
	IDENTICAL=0;

  q1q2=-1;

  nchannels=3;

  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0; // S I=1/2
  ell[1]=0; // S I=3/2
  ell[2]=1; // P I=1/2 (kstar)
  InitWaves();


  nwells=new int[nchannels];
  nwells[0]=2;
  nwells[1]=1;
  nwells[2]=3;

  SquareWell_MakeArrays();

  a[0][0]=0.204619; a[0][1]=0.965437;//
  a[1][0]=0.54948; // 
  a[2][0]=0.23947;  a[2][1]=0.765751; a[2][2]=1.23877;

  V0[0][0]=-9266.34; V0[0][1]=87.0632;//
  V0[1][0]=1057.62; //
  V0[2][0]=-31505.8; V0[2][1]=1391.0; V0[2][2]=-190.325;

  SquareWell_Init();

}

CWaveFunction_kpluspiminus_sqwell::~CWaveFunction_kpluspiminus_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_kpluspiminus_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared;
  double q=GetQ(iq);
  complex<double> psi,psia;
  psia=planewave[iq]->planewave(r,ctheta);
  double x;
  complex<double> DelPhi[3];

  SquareWell_GetDelPhi(iq,r,DelPhi);

  x=q*r/HBARC;

  // For I=1/2 (s & p wave)
  psi=psia;
  psi+=DelPhi[0]/x;
  psi+=3.0*ci*DelPhi[2]/x;
  psisquared=(2.0/3.0)*real(psi*conj(psi));

  // For I=3/2 (s wave only)
  psi=psia;
  psi+=DelPhi[1]/x;
  psisquared+=(1.0/3.0)*real(psi*conj(psi));
  
	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
}

#endif
