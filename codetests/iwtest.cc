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
using namespace std;
#include "constants.h"
#include "wavefunction.h"
#include "misc.h"

int main(){
  const double ROOT2=sqrt(2.0);
  double delta=0.23456;
  int ell=1,ellmax=1,q1q2;
  complex<double> psi,hstar,h,hstarq0,hq0,ci(0.0,1.0);
  double r,x,q=34.0,epsilon=1.0,eta,e1,e2;
  double r0=1.0,rf=4.0,delr=0.001,iw=0.0,iwq0=0.0;
  CPartWave *partwave,*partwaveq0;
  CWaveFunction *wf;
  wf=new CWaveFunction();
  printf("Enter q and delta and ell\n");
  scanf("%lf %lf %d",&q,&delta,&ell);
  //q=35.0;
  //delta=0.4;
  //ell=0;
  q1q2=-1;
  e1=e2=sqrt(q*q+139.57*139.57);
  eta=double(q1q2)*e1*e2/((e1+e2)*137.036*q);
  //eta=1.0E-8;
  printf("eta=%g\n",eta);
  partwave=new CPartWave(eta,q1q2,q,ell,epsilon);
  partwaveq0=new CPartWave(0,0,q,ell,epsilon);
  for(r=r0+0.5*delr;r<rf;r+=delr){
    x=q*r/HBARC;
    hstar=partwave->GetPhiIncoming(r)/x;
    h=partwave->GetPhiOutgoing(r)/x;
    hstarq0=partwaveq0->GetPhiIncoming(r)/x;
    hq0=partwaveq0->GetPhiOutgoing(r)/x;
    /*
      printf("h=(%g,%g), hq0=(%g,%g), e^(ix)/x=(%g,%g)\n",
      real(h),imag(h),real(hq0),imag(hq0),
      real(-ci*ceiphi(x)/x),imag(-ci*ceiphi(x)/x));
      printf("hstar=(%g,%g), hstarq0=(%g,%g)\n",
      real(hstar),imag(hstar),real(hstarq0),imag(hstarq0));
    */
    psi=hstar*Misc::ceiphi(-2.0*delta)+h;
    iw+=0.5*x*x*real(psi*conj(psi))*delr/HBARC;
    psi=hstarq0*Misc::ceiphi(-2.0*delta)+hq0;
    iwq0+=0.5*x*x*real(psi*conj(psi))*delr/HBARC;
  }
  printf("IW=%g =? %g, IWq0=%g =? %g\n",iw,
	 wf->GetIW(ell,r0,q,q1q2,eta,delta)-wf->GetIW(ell,rf,q,q1q2,eta,delta),
	 iwq0,
	 wf->GetIW(ell,r0,q,0,0.0,delta)-wf->GetIW(ell,rf,q,0,0.0,delta));
  delete partwave;

  return 0;
  
}
