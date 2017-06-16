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
#ifndef __INCLUDE_PIPLUSPIPLUS_SQWELL_CC__
#define __INCLUDE_PIPLUSPIPLUS_SQWELL_CC__
#include "wavefunction.h"

using namespace std;

CWaveFunction_pipluspiplus_sqwell::CWaveFunction_pipluspiplus_sqwell(string  parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);
  m1=MPI; 
  m2=MPI;
	IDENTICAL=1;

  q1q2=1;

  nchannels=2;

  ellmax=2;
  InitArrays();
  printf("Arrays Initialized\n");

  ell[0]=0; // S2
  ell[1]=2; // D2

  InitWaves();
  printf("Partial Waves Initialized\n");
  nwells=new int[nchannels];
  nwells[0]=1;
  nwells[1]=1;

  SquareWell_MakeArrays();

  a[0][0]=0.188947;
  a[1][0]=0.452495;
	
  V0[0][0]=52894.6;
  V0[1][0]=16468.6;
  
  SquareWell_Init();

}

CWaveFunction_pipluspiplus_sqwell::~CWaveFunction_pipluspiplus_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_pipluspiplus_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  const double ROOT2=sqrt(2.0);
  double psisquared;
  double q=GetQ(iq);
  complex<double> psi,psia, psib, psisymm, psianti;
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(psia+psib)/ROOT2;
  psianti=(psia-psib)/ROOT2;
  double x;
  complex<double> DelPhi[5];

  SquareWell_GetDelPhi(iq,r,DelPhi);

  x=q*r/HBARC;
  // For I=2
  psi=psisymm;
  psi+=ROOT2*DelPhi[0]/x;
  psi-=ROOT2*5.0*SpherHarmonics::legendre(2,ctheta)*DelPhi[1]/x;
  psisquared=real(psi*conj(psi));

	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
}

#endif
