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
#include "wavefunction.h"

int main(){
//  const double PI=3.14159265358979323844;
  const double HBARC=197.3269602,ROOT2=sqrt(2.0);
  int ell=0,ellmax=10,q1q2;
  complex<double> psi_partwave(0.0,0.0),psi_planewave(0.0,0.0),ci(0.0,1.0);
  complex<double> cipow,phi_in,phi_out,phi1,phi2;
  double sigma;
  double r,x,q=34.0,epsilon=0.1,eta,e1,e2,ctheta;
  double r0=0.2,rf=1,delr=0.1;
  CPartWave **partwave;
  CPlaneWave *planewave;
  //printf("Enter q and q1q2\n");
  //scanf("%lf %d",&q,&q1q2);
  q=350.0; q1q2=1;
  e1=e2=sqrt(q*q+139.57*139.57);
  //eta=double(q1q2)*e1*e2/((e1+e2)*137.036*q);
  eta=0.1;
  printf("eta=%g\n",eta);
  partwave=new CPartWave *[ellmax+1];
  for(ell=0;ell<=ellmax;ell++)
    partwave[ell]=new CPartWave(eta,q1q2,q,ell,epsilon);
  planewave=new CPlaneWave(eta,q1q2,q);
  

  for(r=r0+0.5*delr;r<rf;r+=delr){
    ctheta=0.3456;
    x=q*r/HBARC;
    psi_planewave=planewave->planewave(r,ctheta);
    psi_partwave=0.0;
    cipow=-ci;
    for(ell=0;ell<=ellmax;ell++){
      cipow*=ci;
      sigma=partwave[ell]->sigma;
      //printf("sigma(ell=%d)=%g, ",ell,sigma);
      phi_in=CoulWave::CWincoming(ell,x,eta);
      phi_out=CoulWave::CWoutgoing(ell,x,eta);
      psi_partwave+=0.5*(2.0*ell+1.0)
	*cipow*(partwave[ell]->GetPhiIncoming(r)+
		partwave[ell]->GetPhiOutgoing(r))
	*SpherHarmonics::legendre(ell,ctheta)/x;
    }	   
    //printf("\n");

    printf("____r=%g, x=%g, e^(iqdotr)=(%g,%g) ___\n psi_planewave=(%g,%g), psi_partwave=(%g,%g)\n",
	   r,x,cos(x*ctheta),sin(x*ctheta),
	   real(psi_planewave),imag(psi_planewave),
	   real(psi_partwave),imag(psi_partwave));
  }
  delete planewave;
  for(ell=0;ell<=ellmax; ell++) delete partwave[ell];
  delete [] partwave;
  
  return 0;
  
}
