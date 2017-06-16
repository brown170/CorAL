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
#ifndef __INCLUDE_KPLUSPIPLUS_CC__
#define __INCLUDE_KPLUSPIPLUS_CC__
#include "wavefunction.h"


using namespace std;

CWaveFunction_kpluspiplus_sqwell::CWaveFunction_kpluspiplus_sqwell(string parsfilename) : CWaveFunction(){
  ParsInit(parsfilename);

  m1=MPI; 
  m2=MKAON;
	IDENTICAL=0;

  q1q2=1;

  nchannels=1;

  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");
  ell[0]=0; // S I=3/2
  InitWaves();

  nwells=new int[nchannels];
  nwells[0]=1;

  SquareWell_MakeArrays();

  a[0][0]=1.26401;
  V0[0][0]=49.8849;
  
  SquareWell_Init();

}

CWaveFunction_kpluspiplus_sqwell::~CWaveFunction_kpluspiplus_sqwell(){
  SquareWell_DeleteArrays();
}


double CWaveFunction_kpluspiplus_sqwell::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared;
  double q=GetQ(iq);
  complex<double> psi;
  psi=planewave[iq]->planewave(r,ctheta);
  double x;
  complex<double> DelPhi[1];

  SquareWell_GetDelPhi(iq,r,DelPhi);

  x=q*r/HBARC;
  // For I=3/2
  psi+=DelPhi[0]/x;
  psisquared=real(psi*conj(psi));
	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
}

#endif
