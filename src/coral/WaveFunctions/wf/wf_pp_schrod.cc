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
#ifndef __INCLUDE_PP_SCHROD_CC__
#define __INCLUDE_PP_SCHROD_CC__
#include "wavefunction.h"
using namespace std;

CWaveFunction_pp_schrod::CWaveFunction_pp_schrod(string parsfilename) : CWaveFunction(){
  int ichannel,iq;
  ParsInit(parsfilename);
  m1=MPROTON; 
  m2=m1; 
	IDENTICAL=1;
  q1q2=1;
	
  nchannels=4;
  ellmax=1;
  InitArrays();
  printf("Arrays Initialized\n");
  ell[0]=0;
  ell[1]=ell[2]=ell[3]=1;
	
  InitWaves();
  printf("Partial Waves Initialized\n");
	
  rmax_schrod=6.0;
  nrmax_schrod=600;
  delpsi=new complex<double> **[nchannels];
  for(ichannel=0;ichannel<nchannels;ichannel++){
    delpsi[ichannel]=new complex<double> *[nqmax];
    for(iq=0;iq<nqmax;iq++){
      delpsi[ichannel][iq]=new complex<double>[nrmax_schrod+1];
      schrodinger(ichannel,iq);
    }
  }
}


CWaveFunction_pp_schrod::~CWaveFunction_pp_schrod(){
  int ichannel,iq;
  for(ichannel=0;ichannel<nchannels;ichannel++){
    for(iq=0;iq<nqmax;iq++){
      delete [] delpsi[ichannel][iq];
    }
    delete [] delpsi[ichannel];
  }
  delete [] delpsi;
}

double CWaveFunction_pp_schrod::CalcPsiSquared(int iq,double r,double ctheta){
  int ir1,ir2;
  double w1,w2,psisquared,theta,x,delr,q;
  const double ROOT2=sqrt(2.0);
  complex<double> delpsiS0=0.0,delpsiS1J0=0.0,delpsiS1J1=0.0,delpsiS1J2=0.0;
  complex<double> psia,psib,psisymm,psianti,Xlm00,Xlm10,Xlm11;
  complex<double> psi_samespin,psi_diffspin;
  q=(0.5+iq)*delq;
  delr=rmax_schrod/double(nrmax_schrod);
  if(r<rmax_schrod){
    ir1=lrint(floor(r/delr));
    ir2=ir1+1;
    w2=(r-ir1*delr)/delr;
    w1=1.0-w2;
    delpsiS0=w1*delpsi[0][iq][ir1]+w2*delpsi[0][iq][ir2];
    delpsiS1J0=w1*delpsi[1][iq][ir1]+w2*delpsi[1][iq][ir2];
    delpsiS1J1=w1*delpsi[2][iq][ir1]+w2*delpsi[2][iq][ir2];
    delpsiS1J2=w1*delpsi[3][iq][ir1]+w2*delpsi[3][iq][ir2];
  }
  else{
    delpsiS0=(exp(-2.0*ci*delta[0][iq])-1.0)*partwave[0][iq]->GetPhiIncoming(r);
    delpsiS1J0=(exp(-2.0*ci*delta[1][iq])-1.0)*partwave[1][iq]->GetPhiIncoming(r);
    delpsiS1J1=(exp(-2.0*ci*delta[2][iq])-1.0)*partwave[1][iq]->GetPhiIncoming(r);
    delpsiS1J2=(exp(-2.0*ci*delta[3][iq])-1.0)*partwave[1][iq]->GetPhiIncoming(r);
  }
  theta=acos(ctheta);
  x=q*r/HBARC;
	
  //Notation: Xlm_{l,m_l}
  Xlm00=0.5*ROOT2*sqrt(4.0*PI)*SpherHarmonics::Ylm(0,0,theta,0.0)/x;
  Xlm10=0.5*ROOT2*ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,0,theta,0.0)/x;
  Xlm11=0.5*ROOT2*ci*sqrt(12.0*PI)*SpherHarmonics::Ylm(1,1,theta,0.0)/x;
  psia=planewave[iq]->planewave(r,ctheta);
  psib=planewave[iq]->planewave(r,-ctheta);
  psisymm=(1.0/ROOT2)*(psia+psib);
  psianti=(1.0/ROOT2)*(psia-psib);
	
  psisquared=0.0;
  //S=0
  psi_samespin=psisymm;
  psi_samespin+=delpsiS0*Xlm00;
  psisquared+=0.25*real(psi_samespin*conj(psi_samespin));
  //S=1,m_s=1
  psi_samespin=psianti;
  psi_samespin+=(0.5*delpsiS1J2+0.5*delpsiS1J1)*Xlm10;
  psi_diffspin=(0.5*delpsiS1J2-0.5*delpsiS1J1)*Xlm11;
  psisquared+=0.5*real(psi_samespin*conj(psi_samespin)
		+psi_diffspin*conj(psi_diffspin));
  //S=1,m_s=0
  psi_samespin=psianti;
  psi_samespin+=((2.0/3.0)*delpsiS1J2+(1.0/3.0)*delpsiS1J0)*Xlm10;
  psi_diffspin=((1.0/3.0)*delpsiS1J2-(1.0/3.0)*delpsiS1J0)*Xlm11;
  psisquared+=0.25*real(psi_samespin*conj(psi_samespin))
	+0.5*real(psi_diffspin*conj(psi_diffspin));
	
	psisquared*=RelativisticCorrection(r,iq);
  return psisquared;
}

void CWaveFunction_pp_schrod::schrodinger(int ichannel,int iq){
  char type[10];
  char pname[10];
	
  double q=(iq+0.5)*delq,q2,r;
  double delr=0.001;
  double x1,mu,E,delx2,r0,r1,r2,vr,v12,v22,v11,phase,phase0,sigma;
  double phaseb,phaseb0;
  complex<double> cg;
  int nr=lrint(rmax_schrod/delr);
  int ir,J,S,L=ell[ichannel];
  complex<double> *psiout,*psiout0;
  sprintf(type,"PP");
  if(ichannel==0){
    J=0; S=0; L=0;
    sprintf(pname,"1S0");
  }
  else{
    S=1; L=1; 
    if(ichannel==1){
      J=0;
      sprintf(pname,"3P0");
    }
    else if(ichannel==2){
      J=1;
      sprintf(pname,"3P1");	
    }
    else{
      J=2;
      sprintf(pname,"3C2");
    }
  }	
  E=2.0*sqrt(q*q+m1*m1);
  mu=0.25*E;
  q2=q*q;
  delx2=q2*delr*delr/(HBARC*HBARC);
  cg=CoulWave::cgamma(L+1.0+ci*eta[iq]);
  sigma=atan2(imag(cg),real(cg));
  psiout=new complex<double> [nr+1];
  psiout0=new complex<double> [nr+1];
  complex<double>psi2,psi1,psi0;
	
  r2=nr*delr;
  r1=(nr-1)*delr;
  r0=(nr-2)*delr;
  psiout[nr]=CoulWave::CWoutgoing(L,r2*q/HBARC,eta[iq])*exp(-ci*sigma);
  psiout[nr-1]=CoulWave::CWoutgoing(L,r1*q/HBARC,eta[iq])*exp(-ci*sigma);
  psiout0[nr]=psiout[nr];
  psiout0[nr-1]=psiout[nr-1];
	
  for(ir=nr-2;ir>=1;ir--){
    r1=(ir+1)*delr;
    WaveFunctionRoutines::creid93(&r1,pname,type,&v11,&v12,&v22);
    //v11=Vreid(r1,ichannel);
    vr=2.0*mu*v11/q2;
    x1=r1*q/HBARC;
    psiout0[ir]=2.0*psiout0[ir+1]-psiout0[ir+2]
		+delx2*(-1.0+(2.0*eta[iq]/x1)+double(L)*double(L+1)/(x1*x1))
		*psiout0[ir+1];
    psiout[ir]=2.0*psiout[ir+1]-psiout[ir+2]
		+delx2*(-1.0+vr+(2.0*eta[iq]/x1)+double(L)*double(L+1)/(x1*x1))
		*psiout[ir+1];
    /*if(L==0 && ir<10){
      printf("r=%g, x=%g, vr=%g, vc=%g, psiout0=(%g,%g), psiout=(%g,%g)\n",
				r1-delr,x1-delr*q/HBARC,vr,2.0*eta[iq]/x1,
				real(psiout0[ir]),imag(psiout0[ir]),real(psiout[ir]),imag(psiout[ir]));	
		}*/
  }
  phase=atan2(imag(psiout[1]),real(psiout[1]))+0.5*PI;
  phase0=atan2(imag(psiout0[1]),real(psiout0[1]))+0.5*PI;
  phaseb=atan2(imag(psiout[2]),real(psiout[2]))+0.5*PI;
  phaseb0=atan2(imag(psiout0[2]),real(psiout0[2]))+0.5*PI;
  phase=2.0*phase-phaseb;
  phase0=2.0*phase0-phaseb0;
  delta[ichannel][iq]=-phase+phase0;
  //printf("channel %d: L=%d,  q=%g, phase=%g, phase0=%g, sigma=%g, delta=%g\n",
		// ichannel,L, q,phase,phase0,sigma,180.0*delta[ichannel][iq]/PI);
	
  int ir1,ir2,ir_schrod;
  double w1,w2;
  delpsi[ichannel][iq][0]=0.0;
  for(ir_schrod=1;ir_schrod<=nrmax_schrod;ir_schrod++){
    r=rmax_schrod*ir_schrod/double(nrmax_schrod);
    ir1=lrint(floor(r/delr));
    ir2=ir1+1;
    if(ir2<=nr){
      w1=(ir2*delr-r)/delr;
      w2=1.0-w1;
      delpsi[ichannel][iq][ir_schrod]
			=w2*((psiout[ir2]-psiout0[ir2])
				+conj(exp(-2.0*ci*phase)*psiout[ir2]
					-exp(-2.0*ci*phase0)*psiout0[ir2]));
      delpsi[ichannel][iq][ir_schrod]
			+=w1*((psiout[ir1]-psiout0[ir1])
	      +conj(exp(-2.0*ci*phase)*psiout[ir1]
					-exp(-2.0*ci*phase0)*psiout0[ir1]));
    }
    else if(ir1==nr)
      delpsi[ichannel][iq][ir_schrod]
			=((psiout[ir1]-psiout0[ir1])
				+conj(exp(-2.0*ci*phase)*psiout[ir1]
					-exp(-2.0*ci*phase0)*psiout0[ir1]));
			else{
				printf("OOPS, out of bounds, shouldn't be here, ir1=%d, nr=%d\n",ir1,nr);
				exit(1);
			}
			/*
      if(iq==0 && ichannel==0){
				x=q*r/HBARC;
				printf("r=%g, ir1=%d, ir2=%d, w1=%g, w2=%g, delpsi=(%g,%g)\n",r,ir1,ir2,w1,w2,
					real(delpsi[ichannel][iq][ir_schrod]),
					imag(delpsi[ichannel][iq][ir_schrod]));
      }*/
  }
	
  delete [] psiout;
  delete [] psiout0;
	
}

double CWaveFunction_pp_schrod::Vreid(double r,int ipart){
  double pmux,f1,f2,f4,f6,f7,vr=0.0;
  /* See the appendix of B.D. Day, PRC24, p. 1203 (1981).
	with Tensor forces neglected */
  if(ipart==0){
    /* S=0 */
    pmux=r*0.7;
    f1=exp(-pmux);
    f4=(f1*f1*f1*f1);
    f7=f4*(f1*f1*f1);
    vr=-10.463*f1/pmux-1650.6*f4/pmux+6484.2*f7/pmux;
  }
  if(ipart>0){
    /* S=1 */
    pmux=r*0.7;
    f1=exp(-pmux);
    f2=f1*f1;
    f4=f2*f2;
    f6=f4*f2;
    vr=((-10.463/3.0)*f1-933.48*f4+4152.1*f6)/pmux;
  }
  return vr;
}

#endif
