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
#ifndef __CWAVEFUNCTION_WF_PIPI_PHASESHIFT_CC__
#define __CWAVEFUNCTION_WF_PIPI_PHASESHIFT_CC__
#include "wavefunction.h"

void WaveFunctionRoutines::getphaseshift_pipi(int I,int ell,double q,double *delta,double *ddeltadq){
  double a,a2,a3,reff,x,R,tandelta;
  const double MPI=139.58,MRHO=771.1,Gamma=149.2; // from pdg
  const double qrho=sqrt(0.25*MRHO*MRHO-MPI*MPI);
  double q1,q2,delta1,delta2,M,piself,delq;
  double ff,mscreen;

  if((ell+I)%2 != 0){
    printf("pipi phase shift: bad (ell,I) comb, ell=%d, I=%d\n",ell,I);
    exit(1);
  }

  // For I=0, See R. Kaminski et al. Acta Phys. Polonica, B31 (2000)
  // M.E. Sevoir, NPA 543, 275c (1992).
  // G. Grayer et al NPB75, 189 (1975)
  // L. Rosselet, et al, PRD, 574 (1977).
  // P. Estabrooks and A.D. Martin, NPB79, p. 301 (1974)
  // The I=0,ell=0 formula will only work for energies below the 2K threshold
  if(I==0){
    if(ell==0){
      // This uses the scattering length of 0.204/MPI for the q->0 solution,
      // and adds an extra piece proportional to (M-2MPI) to 
      // roughly fit the higher E data.
      a=0.204/MPI;
      a2=290.0;
      a3=625.0;
      tandelta=a*q/(1.0-(q/a2)+pow(q/a3,2));
      *delta=atan(tandelta);
      *delta=*delta+PI*(fabs(*delta)-*delta)/fabs(2**delta);
      *ddeltadq=pow(cos(*delta),2)*(tandelta/q)*(1.0-q*q/(a3*a3))
	/(1.0-(q/a2)+q*q/(a3*a3));
    }
    if(ell==2){
      // The data for ell=2 aren't too smooth (see Estabrooks).
      // This is an interpolation
      // assuming delta ~ q^5 with delta(MPIpi=1.0 GeV, q=480) = 9 degrees
      a2=1000.0;
      a=9.0*(PI/180.0)/pow(480.0,5);
      x=q/(1.0+q*q/(a2*a2));
      tandelta=a*pow(x,5);
      *delta=atan(tandelta);
      *ddeltadq=pow(cos(*delta),2)*(tandelta/q)
	*((2*ell+1)-(2*ell-1)*q*q/(a2*a2))
	/(1+q*q/(a2*a2));
    }
    if(ell>2){
      *delta=0.0;
      *ddeltadq=0.0;
    }
  }

  // For I=1, phase shifts are dominated by the rho
  // See data from S.D. Protopopescu et al PRD 7, p. 1279 (1973).
  //
  if(I==1){
    if(ell==1){
      mscreen=1200;
      delq=1.0;
      if(q<20.0) delq=q/20.0;
      if(q<0.5*delq) q=0.5*delq;
	
      q1=q-0.5*delq;
      M=2.0*sqrt(q1*q1+MPI*MPI);
      piself=MRHO*Gamma*pow(q1/qrho,3)*(MRHO/M);
      //ff=mscreen*mscreen/(mscreen*mscreen+q1*q1);
      ff=1.0;
      piself=piself*ff;
      delta1=atan2(piself,MRHO*MRHO-M*M);
	
      q2=q+0.5*delq;
      M=2.0*sqrt(q2*q2+MPI*MPI);
      piself=MRHO*Gamma*pow(q2/qrho,3)*(MRHO/M);
      //ff=mscreen*mscreen/(mscreen*mscreen+q2*q2);
      ff=1.0;
      piself=piself*ff;
      delta2=atan2(piself,MRHO*MRHO-M*M);
	
      *delta=0.5*(delta1+delta2);
      *ddeltadq=(delta2-delta1)/delq;
    }
    if(ell>1){
      *delta=0.0;
      *ddeltadq=0.0;
    }
  }

  // For I=2, See M.J. Losty, NPB 69, p. 185 (1974)
  if(I==2){
    if(ell==0){
      a=-0.13;
      reff=1.0;
      /* cotdelta=(HBARC/(q*a))+0.5*q*reff/HBARC;
       *delta=atan(1.0/cotdelta);
       *ddeltadq=pow(sin(*delta),2)*(-0.5*reff/HBARC+HBARC/(q*q*a)); */
      
      R=sqrt(0.5*fabs(a*reff));
      x=q/(1.0+pow(q*R/HBARC,2));
      tandelta=(a/HBARC)*x;
      *delta=atan(tandelta);
      *ddeltadq=(pow(cos(*delta),2)*tandelta/q)*(2.0*ell+1)
	*(1.0-pow(q*R/HBARC,2))/(1.0+pow(q*R/HBARC,2));
    }
    if(ell==2){
      /* *delta=-8.4*pow(q/1000.0,5)+12.5*pow(q/1000.0,6);
       *ddeltadq=-8.4*(5.0/q)*pow(q/1000.0,3)+12.5*(6.0/q)*pow(q/1000.0,4); */
      R=0.5;
      x=q/(1.0+pow(q*R/HBARC,2));
      tandelta=-(8.4/pow(1000.0,5))*pow(x,5);
      *delta=atan(tandelta);
      *ddeltadq=(pow(cos(*delta),2)*tandelta/q)*(2.0*ell+1)
	*(1.0-pow(q*R/HBARC,2))/(1.0+pow(q*R/HBARC,2));
    }
    if(ell>2){
      *delta=0.0;
      *ddeltadq=0.0;
    }
  }

}

#endif
