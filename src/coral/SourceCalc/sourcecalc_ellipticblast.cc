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
using namespace std;

#include "sourcecalc.h"

CSourceCalc_EllipticBlast::CSourceCalc_EllipticBlast(){
  InitSPars();
  randy=new CRandom(-1234);
}

void CSourceCalc_EllipticBlast::InitSPars(){
  // DEFAULT VALUES
  parameter::set(spars,"Rx",13);
  parameter::set(spars,"Ry",13);
  parameter::set(spars,"Tau",12);
  parameter::set(spars,"BetaX",0.7);
  parameter::set(spars,"BetaY",0.7);
  parameter::set(spars,"T",110);
  parameter::set(spars,"Pt",600);
  parameter::set(spars,"Phi",0.0);
  parameter::set(spars,"EtaG",2.0);  
  parameter::set(spars,"Ma",938.28);
  parameter::set(spars,"Mb",139.58);
  parameter::set(spars,"Nsample",1000);
}

void CSourceCalc_EllipticBlast::SetSPars(double Rxset,double Ryset,double Tauset,
				double BetaXset,double BetaYset,
				double Tset,double Ptset,
				double Phiset,double EtaGset,
				double Maset,double Mbset){
  InitSPars();
  parameter::set(spars,"Rx",Rxset);
  parameter::set(spars,"Ry",Ryset);
  parameter::set(spars,"Tau",Tauset);
  parameter::set(spars,"BetaX",BetaXset);
  parameter::set(spars,"BetaY",BetaYset);
  parameter::set(spars,"T",Tset);
  parameter::set(spars,"Pt",Ptset);
  parameter::set(spars,"Phi",Phiset);
  parameter::set(spars,"EtaG",EtaGset);  
  parameter::set(spars,"Ma",Maset);  
  parameter::set(spars,"Mb",Mbset);

}

void CSourceCalc_EllipticBlast::SetSPars(double Rset, double Tauset, double Betaset, double Tset, double Ptset){

  InitSPars();
  parameter::set(spars,"Rx",Rset);
  parameter::set(spars,"Ry",Rset);
  parameter::set(spars,"Tau",Tauset);
  parameter::set(spars,"BetaX",Betaset);
  parameter::set(spars,"BetaY",Betaset);
  parameter::set(spars,"T",Tset);
  parameter::set(spars,"Pt",Ptset);

}

void CSourceCalc_EllipticBlast::CalcS(CCHArray *A){
  int alpha;
  double pa[4],pb[4],ptot[4];
  double cphi,sphi,ma,mb,phi,Pt,T,tau,volume,delr;
  int nsample;
  double rcm[4],**ra,**rb;
  int ir,imc,ia,ib,nbmax,nrmax;
  double x,y,z,xbar,ybar,zbar,x2bar,y2bar,z2bar,ex,ey,ez,snorm,r;
  bool SameMass;
  const double PI=4.0*atan(1.0);
  delr=A->GetRADSTEP();
  nrmax=A->GetNRADIAL();

  nsample=parameter::getI(spars,"Nsample",1000);
  ma=parameter::getD(spars,"Ma",-999);
  mb=parameter::getD(spars,"Mb",-999);
  SameMass=0;
  if(fabs(ma-mb)<1.0) SameMass=1;
  pa[3]=pb[3]=0.0;
  phi=parameter::getD(spars,"Phi",-999);
  Pt=parameter::getD(spars,"Pt",-999);
  T=parameter::getD(spars,"T",-999);
  tau=parameter::getD(spars,"Tau",-999);
  sphi=sin(phi);
  cphi=cos(phi);
  pa[1]=cphi*Pt*ma/(ma+mb);
  pb[1]=cphi*Pt*mb/(ma+mb);
  pa[2]=sphi*Pt*ma/(ma+mb);
  pb[2]=sphi*Pt*mb/(ma+mb);
  pa[0]=sqrt(ma*ma+pa[1]*pa[1]+pa[2]*pa[2]+pa[3]*pa[3]);
  pb[0]=sqrt(mb*mb+pb[1]*pb[1]+pb[2]*pb[2]+pb[3]*pb[3]);
  for(alpha=0;alpha<4;alpha++) ptot[alpha]=pa[alpha]+pb[alpha];

  ra=new double *[nsample];
  for(imc=0;imc<nsample;imc++) ra[imc]=new double[4];
  Get_r(pa,nsample,ra);
  if(SameMass){
    rb=ra;
  }
  else{
    rb=new double *[nsample];
    for(imc=0;imc<nsample;imc++) rb[imc]=new double[4];
    Get_r(pb,nsample,rb);
  }

  rcm[0]=0.0;
  xbar=ybar=zbar=x2bar=y2bar=z2bar=0.0;
  for(ia=0;ia<nsample;ia++){
    nbmax=nsample;
    if(SameMass) nbmax=ia-1;
    for(ib=0;ib<nbmax;ib++){
      rcm[1]=ra[ia][1]-rb[ib][1];
      rcm[2]=ra[ia][2]-rb[ib][2];
      rcm[3]=ra[ia][3]-rb[ib][3];

      r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
      x=rcm[1];
      y=rcm[2];
      z=rcm[3];
      xbar+=x;
      x2bar+=x*x;
      ybar+=y;
      y2bar+=y*y;
      zbar+=z;
      z2bar+=z*z;

      ir=int(floor(r/delr));
      if(ir<nrmax){
	ex=x/r; ey=y/r; ez=z/r;
	A->IncrementAExpArrayFromE(ex,ey,ez,1.0,ir);
      }
    }
    if(10*(ia+1)%nsample==0)
      printf("finished %g percent\n",100*double(ia+1)/double(nsample));
  }
  A->FillRemainderX();
  if(SameMass) snorm=2.0/double(nsample*(nsample-1));
  else snorm=1.0/double(nsample*nsample);
  xbar*=snorm;
  x2bar*=snorm;
  ybar*=snorm;
  y2bar*=snorm;
  zbar*=snorm;
  z2bar*=snorm;
  x2bar=x2bar-xbar*xbar;
  y2bar=y2bar-ybar*ybar;
  z2bar=z2bar-zbar*zbar;

  printf("xbar=%g, ybar=%g, zbar=%g\n",xbar,ybar,zbar);
  printf("Effective Gaussian Radii: Rout=%g, Rside=%g, Rlong=%g\n",
	 sqrt(0.5*x2bar),sqrt(0.5*y2bar),sqrt(0.5*z2bar));

  for(ir=0;ir<nrmax;ir++){
    volume=(4.0*PI/3)*(pow((ir+1)*delr,3)
		       -pow(double(ir)*delr,3));
    A->ScaleArray(snorm/volume,ir);
  }

  for(imc=0;imc<nsample;imc++){
    delete ra[imc];
    if(!SameMass) delete rb[imc];
  }
  delete ra;
  delete rb;

}

void CSourceCalc_EllipticBlast::Get_r(double *p, int nsample, double **r){
  double phi,eta,arg,etaG,u[4],x,y,z,t,Rx,Ry,weight;
  double uxmax,uymax,betaxmax,betaymax,eprime,eu,T,tau;
  double pt,gamma,gammav,rap,rout,rlong,rside,sinhy,coshy;
  int imc;
  const double PI=4.0*atan(1.0);
  double m=sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
  etaG=parameter::getD(spars,"EtaG",-999);
  Rx=parameter::getD(spars,"Rx",-999);
  Ry=parameter::getD(spars,"Ry",-999);
  T=parameter::getD(spars,"T",-999);
  tau=parameter::getD(spars,"Tau",-999);
  betaxmax=parameter::getD(spars,"BetaX",-999);
  betaymax=parameter::getD(spars,"BetaY",-999);
  uxmax=betaxmax/sqrt(1.0-betaxmax*betaxmax);
  uymax=betaymax/sqrt(1.0-betaymax*betaymax);
  rap=atanh(p[3]/p[0]);
  pt=sqrt(p[1]*p[1]+p[2]*p[2]);
  gammav=pt/m;
  gamma=sqrt(1.0+gammav*gammav);

  for(imc=0;imc<nsample;imc++){
    do{
      eta=etaG*randy->gauss();
    TRY_AGAIN:
      x=(1.0-2.0*randy->ran());
      y=(1.0-2.0*randy->ran());
      if(x*x+y*y>1.0) goto TRY_AGAIN;
      u[1]=uxmax*x;
      u[2]=uymax*y;
      x=x*Rx;
      y=y*Ry;
      u[3]=sinh(eta);
      u[0]=sqrt(1.0+u[1]*u[1]+u[2]*u[2]+u[3]*u[3]);
      
      eprime=p[0]*cosh(eta)-p[3]*u[3];
      eu=u[0]*p[0]-u[1]*p[1]-u[2]*p[2]-u[3]*p[3];

      weight=(eprime/p[0])*exp(-(eu-m)/T);
      if(weight>1.0){
	printf("DISASTER! weight > 1.0\n");
	exit(1);
      }
    } while(weight<randy->ran());

    z=tau*sinh(eta);
    t=tau*cosh(eta);

    rout=(p[1]*x+p[2]*y)/pt;
    rside=(p[1]*y-p[2]*x)/pt;
    sinhy=sinh(rap);
    coshy=cosh(rap);
    rlong=coshy*z-sinhy*t;
    t=coshy*t-sinhy*z;
    r[imc][0]=gamma*t-gammav*rout;
    r[imc][1]=gamma*rout-gammav*t;
    r[imc][2]=rside;
    r[imc][3]=rlong;

  }
}
