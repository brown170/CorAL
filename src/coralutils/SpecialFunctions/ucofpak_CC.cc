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
//  This is a c++ rewrite from FORTRAN.  The functions are the complex
// gamma and confluent hypergeometric function
// 
#include <iostream>
#include <cmath>
#include <complex>
#include "ucofpak_CC.h"
#include "constants.h"

using namespace std;

double_complex cdgamma(double_complex Z)
{
  double_complex CDGAMMA;
  double_complex U,V,H,S,FK;
	const double PI=4.0*atan(1.0);
  double G[16],X;
  int L;

  G[0]=41.624436916439068;
  G[1]=-51.224241022374774;
  G[2]=+11.338755813488977;
  G[3]=-0.747732687772388;
  G[4]=+0.008782877493061;
  G[5]=-0.000001899030264;
  G[6]=+0.000000001946335;
  G[7]=-0.000000000199345;
  G[8]=+0.000000000008433;
  G[9]=+0.000000000001486;
  G[10]=-0.000000000000806;
  G[11]=+0.000000000000293;
  G[12]=-0.000000000000102;
  G[13]=+0.000000000000037;
  G[14]=-0.000000000000014;
  G[15]=+0.000000000000006;

  U=Z;
  X=real(U);
  
  V=1.0-U;
  L=1;
  if (X>=0.0)
    {
      V=U+1.0;
      L=2;
    }
  if (X>=1.0)
    {
      V=U;
      L=3;
    }
  
  H=double_complex(1.0,0.0);
  
  S=G[0];
  FK=0;
  for(int K=1;K<16;K++){
    FK=K-1.0;
    H=((V-(FK+1.0))/(V+FK))*H;
    S=S+G[K]*H;
  }
  
  H=V+4.5;
  CDGAMMA=2.506628274631001*exp((V-0.5)*log(H)-H)*S;
  if(L==1){
    CDGAMMA=PI/(sin(PI*U)*CDGAMMA);
    return CDGAMMA;
  }
  
  if(L==2){
    CDGAMMA=CDGAMMA/U;
    return CDGAMMA;
  }
  if(L==3){
    return CDGAMMA;
  }
  cerr<<"Should not get here"<<endl;
  return 0.0;
}


double_complex cdfhg(double_complex A, double_complex B, double_complex Z)
{

  //C  DOUBLE COMPLEX CONFLUENT HYPERGEOMETRIC FUNCTION
  //C  P. DANIELEWICZ   SEP. 4, 1992, OCT. 30, 1995
  //C  ON THE BASIS OF M. ABRAMOWITZ AND I. A. STEGUN,
  //C  HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER, N.Y., 1972
  //   IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  double_complex CDFHG;
  double_complex AA,ZZ,FZZ;
  double_complex DCON,BA,A1;
  double_complex CI(0.0,1.0);
  const int KMAX=150;
  const int KMAX1=20;
  const double TOL=2e-15;
  const double TOLA=2e-15;

  const double PI=4.0*atan(1.0);
  int IA,IBA,N;
  double RZ;

  IA=int(A.real());
  if(IA<=0 && abs(A-(double_complex)(IA))<TOL)
    {
      FZZ=1.0;
      AA=A;
      ZZ=Z;
      goto label100;                          // reg. series
    }
  
  IBA=int((B-A).real());
  if(IBA<=0 && abs(B-A-(double_complex)(IBA))<TOL)    // Kummer transformation
    {
      FZZ=exp(Z);
      AA=B-A;
      ZZ=-Z;
      goto label100;
    }
  
  RZ=real(Z);
  if(RZ>0.0)
    {
      FZZ=1.0;
      AA=A;
      ZZ=Z;
    }else{           // Kummer
      FZZ=exp(Z);
      AA=B-A;
      ZZ=-Z;
    }
  
  if(abs(AA)<TOL)                      // no sense to sum
    {
      CDFHG=1.0;
      goto label200;
    }
  
  if(abs(AA-B)<TOL*abs(AA))         // sum turns into exp
    {
      CDFHG=exp(ZZ);
      goto label200;
    }

  if(abs(ZZ)<24.0)goto label100;      // deciding on asympt expns

  CDFHG=exp(ZZ)*pow(ZZ,(AA-B))*cdgamma(B)/cdgamma(AA);
  N=1;
  BA=B-AA;
  A1=1.0-AA;
  DCON=CDFHG*BA*A1/ZZ;
  CDFHG=CDFHG+DCON;

  for(int K=2;K<KMAX1;K++)
    {
      N=N+1;
      BA=BA+1.0;
      A1=A1+1.0;
      DCON=DCON*BA*A1/((double_complex)(N)*ZZ);
      CDFHG=CDFHG+DCON;
      if(abs(DCON)<=TOLA*abs(CDFHG))goto label40;  // 2'nd pt of asmp
    }

  DCON=DCON*(-.5+(.125+.25*B-.5*AA+.25*(ZZ-(double_complex)(KMAX1)))/ZZ);
  CDFHG=CDFHG+DCON;
label40:   

  DCON=exp(CI*PI*AA)/pow(ZZ,AA)*cdgamma(B)/cdgamma(B-AA);
  CDFHG=CDFHG+DCON;
  N=1;
  BA=1.0+AA-B;
  A1=AA;
  DCON=DCON*A1*BA/(-ZZ);
  CDFHG=CDFHG+DCON;

  for(int K=2;K<KMAX1;K++)
    {
      N=N+1;
      BA=BA+1.0;
      A1=A1+1.0;
      DCON=DCON*BA*A1/(-(double_complex)(N)*ZZ);
      CDFHG=CDFHG+DCON;
      if(abs(DCON)<=TOLA*abs(CDFHG))goto label200;  // exit
    }
  DCON=DCON*(-1.0/3.0-B+AA+AA+ZZ-(double_complex)(KMAX1));
  CDFHG=CDFHG+DCON;
  
  goto label200;

label100:

  N=1;
  BA=B;
  A1=AA;
  DCON=A1*ZZ/BA;
  CDFHG=1.0+DCON;
  for(int K=2;K<KMAX;K++)
    {
      N=N+1;
      BA=BA+1.0;
      A1=A1+1.0;
      DCON=DCON*A1*ZZ/((double_complex)(N)*BA);
      CDFHG=CDFHG+DCON;
      if(abs(DCON)<=TOL*abs(CDFHG))goto label200;
    }

  cerr<<"NO CONVERGENCE IN CDFHG"<<endl;

label200:  
      CDFHG=CDFHG*FZZ;
      return CDFHG;
}

/*


      FUNCTION CDFHGV(N,A,B,Z)
      COMPLEX*16 CDFHGV,CDFHG,A,B,Z,AA,BB,AN,BN
C  N'TH DERIVATIVE OF DOUBLE COMPLEX CONFLUENT HYPERGEOMETRIC FUNCTION
C  P. DANIELEWICZ, NOV. 16, 1995
C  ON THE BASIS OF M. ABRAMOWITZ AND I. A. STEGUN,
C  HANDBOOK OF MATHEMATICAL FUNCTIONS, DOVER, N.Y., 1972
C
      AA=A
      AN=A
      BB=B
      BN=B
      DO I=1,N-1
        AA=AA+1D0
        AN=AN*AA
        BB=BB+1D0
        BN=BN*BB
      ENDDO
      CDFHGV=AN/BN*CDFHG(A+N,B+N,Z)
C
      END
*/
