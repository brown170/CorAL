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
#ifndef __CWAVEFUNCTION_WF_GENERIC_CC__
#define __CWAVEFUNCTION_WF_GENERIC_CC__
#include "wavefunction.h"

CWaveFunction_generic::CWaveFunction_generic(
					     string  parsfilename,
					     int q1q2set,double m1set,
					     double m2set,
					     double symmweightset): CWaveFunction(){
  generic=1;
  ParsInit(parsfilename);
  m1=m1set;
  m2=m2set;
  muscale=m1*m2/(m1+m2);
  mu=muscale;
  symmweight=symmweightset;
  q1q2scale=q1q2set;
  q1q2=q1q2scale;
  nchannels=0;
  ellmax=0;
  InitArrays();
  InitWaves();
  printf("initialization finished\n");
}

void CWaveFunction_generic::reset(int q1q2set,double m1set,double m2set,
				  double symmweightset){
  m1=m1set;
  m2=m2set;
  mu=m1*m2/(m1+m2);
  if(q1q2*q1q2set<0){
    printf("Illegal: Trying to reset q1q2 to opposite charge\n");
    exit(1);
  }
  q1q2=q1q2set;
  symmweight=symmweightset;
}

double CWaveFunction_generic::CalcPsiSquared(int iq,double r,double ctheta){
  double psisquared,asymmweight;
  complex<double> psi1,psi2,psisymm,psiasymm;
  const double ROOT2=sqrt(2.0);

  if(iq>=nqmax){
    printf("iq too large!\n");
    psisquared=1.0;
  }
  else{
    psi1=planewave[iq]->planewave(r,ctheta);
    psi2=conj(psi1);
    
    asymmweight=1.0-symmweight;
    psisymm=(psi1+psi2)/ROOT2;
    psiasymm=(psi1-psi2)/ROOT2;
    
    psisquared=symmweight*real(psisymm*conj(psisymm))
      +asymmweight*real(psiasymm*conj(psiasymm));
  }
	if(q1q2!=0) psisquared*=RelativisticCorrection(r,iq);
  return psisquared;

}

#endif
