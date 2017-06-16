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
#ifndef __INCLUDE_PARTWAVE_CC__
#define __INCLUDE_PARTWAVE_CC__
#include "wavefunction.h"

CPartWave::CPartWave(double etaset,int q1q2set,double qset,int ellset,
		     double epsilonset){
  epsilon=epsilonset;
  delr=0.01;
  nrmax=3000;
  complex<double> cg;
  ci=complex<double>(0.0,1.0);
  q=qset;
  q1q2=q1q2set;
  eta=etaset;
  ell=ellset;
  sigma=0.0;
  if(q1q2!=0){
    cg=CoulWave::cgamma(ell+1.0+ci*eta);
    sigma=atan2(imag(cg),real(cg));
  }
  phi_init();
}

CPartWave::CPartWave(double etaset,int q1q2set,double qset,int ellset,
		     double epsilonset,int nrmaxset,double delrset){
  epsilon=epsilonset;
  nrmax=nrmaxset;
  complex<double> cg;
  ci=complex<double>(0.0,1.0);
  q=qset;
  q1q2=q1q2set;
  eta=etaset;
  ell=ellset;
  sigma=0.0;
  if(q1q2!=0){
    cg=CoulWave::cgamma(ell+1.0+ci*eta);
    sigma=atan2(imag(cg),real(cg));
  }
  phi_init();
}

CPartWave::~CPartWave(){
  delete [] phi;
};

void CPartWave::phi_init(){
  complex<double> phi0,phi1,phi2,phitest;
	int ir;
  double x,r;

  phi=new complex<double> [nrmax+1];

  if(q1q2!=0){
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      x=q*r/HBARC;
      phi[ir]=CoulWave::CWincoming(ell,x,eta)*Misc::ceiphi(-sigma);
    }
    
  }
  else{
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      x=q*r/HBARC;
      phi[ir]=x*Bessel::hstarn(ell,x);
    }
  }

}

complex<double> CPartWave::GetPhiIncoming(double r){
  double x,a;
  complex<double> answer;
  int ir;
  ir=int(floor(r/delr));
  if(ir<nrmax){
    a=(r-delr*(ir+0.5))/delr;
    if(a>0.0 || ir==0){
      answer=(1.0-a)*phi[ir]+a*phi[ir+1];
    }
    else{
      answer=(1.0+a)*phi[ir]-a*phi[ir-1];
    }
  }
  else{
    x=q*r/HBARC;
    if(q1q2!=0){
      answer=CoulWave::CWincoming(ell,x,eta);
      answer=answer*Misc::ceiphi(-sigma);
    }
    else
      answer=x*Bessel::hstarn(ell,x);
  }
  return answer;
}

complex<double> CPartWave::GetPhiOutgoing(double r){
  complex<double> answer;
  answer=conj(GetPhiIncoming(r));
  if(q1q2!=0) answer=answer*Misc::ceiphi(-2.0*sigma);
  return answer;
}

#endif
