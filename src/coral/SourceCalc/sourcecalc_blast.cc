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
#ifndef __INCLUDE_SOURCECALC_BLAST__
#define __INCLUDE_SOURCECALC_BLAST__
#include "sourcecalc.h"

using namespace std;

CSourceCalc_Blast::CSourceCalc_Blast(){
	InitSPars();
	randy=new CRandom(-1234);
}

void CSourceCalc_Blast::InitSPars(){
	// DEFAULT VALUES
	parameter::set(spars,"lambda",0.5);
	parameter::set(spars,"R",13.0);
	parameter::set(spars,"Tau",12.0);
	parameter::set(spars,"DelTau",3.0);
	parameter::set(spars,"Beta",0.7);
	parameter::set(spars,"T",110.0);
	parameter::set(spars,"Pt",600.0);
	parameter::set(spars,"Phi",0.0);
	parameter::set(spars,"EtaG",1.5);  
	parameter::set(spars,"Ma",938.28);
	parameter::set(spars,"Mb",139.58);
	parameter::set(spars,"Nsample",1000);
}

void CSourceCalc_Blast::SetSPars(double lambdaset,
	double Rset,double Tauset,double DelTauset,
	double Betaset,double Tset,double Ptset,
double EtaGset,double Maset,double Mbset){
	InitSPars();
	parameter::set(spars,"lambda",lambdaset);
	parameter::set(spars,"R",Rset);
	parameter::set(spars,"Tau",Tauset);
	parameter::set(spars,"Beta",Betaset);
	parameter::set(spars,"DelTau",DelTauset);
	parameter::set(spars,"T",Tset);
	parameter::set(spars,"Pt",Ptset);
	parameter::set(spars,"EtaG",EtaGset);  
	parameter::set(spars,"Ma",Maset);  
	parameter::set(spars,"Mb",Mbset);

}

void CSourceCalc_Blast::SetSPars(double lambdaset, double Rset, double Tauset, double DelTauset){
	InitSPars();
	parameter::set(spars,"lambda",lambdaset);
	parameter::set(spars,"R",Rset);
	parameter::set(spars,"Tau",Tauset);
	parameter::set(spars,"DelTau",DelTauset);
}

void CSourceCalc_Blast::CalcS(CCHArray *A){
	int nsample;
	double rcm[4],*ra,*rb;
	double delr,ma,mb,volume;
	int ir,ia,ib,nbmax,nrmax;
	double x,y,z,xbar,ybar,zbar,x2bar,y2bar,z2bar,ex,ey,ez,snorm,r,lambda;
	bool SameMass;
	CMCList *lista,*listb;
	delr=A->GetRADSTEP();
	nrmax=A->GetNRADIAL();

	ma=parameter::getD(spars,"Ma",-999);
	mb=parameter::getD(spars,"Mb",-999);
	nsample=parameter::getI(spars,"Nsample",1000);
	lambda=parameter::getD(spars,"lambda",-999);
	SameMass=0;
	if(fabs(ma-mb)<1.0) SameMass=1;

	lista=new CMCList(nsample);
	if(!SameMass) listb=new CMCList(nsample);
	else listb=lista;
	CalcS(lista,listb);

	rcm[0]=0.0;
	xbar=ybar=zbar=x2bar=y2bar=z2bar=0.0;
	for(ia=0;ia<nsample;ia++){
		nbmax=nsample;
		if(SameMass) nbmax=ia-1;
		for(ib=0;ib<nbmax;ib++){
			ra=lista->GetR(ia);
			rb=listb->GetR(ib);

			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
			x=rcm[1];
			y=rcm[2];
			z=rcm[3];
#ifdef __CALC_SCBLAST_GAUSSPARS__
			xbar+=x;
			x2bar+=x*x;
			ybar+=y;
			y2bar+=y*y;
			zbar+=z;
			z2bar+=z*z;
#endif
			ir=int(floor(r/delr));
			if(ir<nrmax){
				ex=x/r; ey=y/r; ez=z/r;
				A->IncrementAExpArrayFromE(ex,ey,ez,lambda,ir);
			}
		}
		//if(10*(ia+1)%nsample==0)
		//printf("finished %g percent\n",100*double(ia+1)/double(nsample));
	}
	A->FillRemainderX();
	if(SameMass) snorm=2.0/double(nsample*(nsample-1));
	else snorm=1.0/double(nsample*nsample);

#ifdef __CALC_SCBLAST_GAUSSPARS__
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
#endif

	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)
			-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}

	delete lista;
	if(!SameMass) delete listb;

}

void CSourceCalc_Blast::CalcS(CMCList *lista, CMCList *listb){
	double pa[4],pb[4];
	double ma,mb,Pt,lambda;
	int nsample;

	lambda=parameter::getD(spars,"lambda",1.0);
	lista->SetNorm(sqrt(lambda));
	listb->SetNorm(sqrt(lambda));
	nsample=parameter::getI(spars,"Nsample",1000);
	ma=parameter::getD(spars,"Ma",-999);
	mb=parameter::getD(spars,"Mb",-999);
	pa[3]=pb[3]=0.0;
	Pt=parameter::getD(spars,"Pt",-999);
	pa[1]=Pt*ma/(ma+mb);
	pb[1]=Pt*mb/(ma+mb);
	pa[2]=pa[3]=0.0;
	pb[2]=pb[3]=0.0;
	pa[0]=sqrt(ma*ma+pa[1]*pa[1]);
	pb[0]=sqrt(mb*mb+pb[1]*pb[1]);

	if(lista->GetNMC()!=nsample) lista->Resize(nsample);
	GetMCList(pa,lista);

	if(lista!=listb){
		if(listb->GetNMC()!=nsample) listb->Resize(nsample);
		GetMCList(pb,listb);
	}

}

void CSourceCalc_Blast::GetMCList(double *p, CMCList *mclist){
	double eta,etaG,u[4],x,y,z,t,tt,R,weight;
	double umax,betamax,eprime,eu,eumin,T,tau,tau0,deltau;
	double pt,gamma,gammav;
	int imc,nsample;
	double m=sqrt(p[0]*p[0]-p[1]*p[1]);
	etaG=parameter::getD(spars,"EtaG",-999);
	R=parameter::getD(spars,"R",-999);
	T=parameter::getD(spars,"T",-999);
	tau0=parameter::getD(spars,"Tau",-999);
	betamax=parameter::getD(spars,"Beta",-999);
	deltau=parameter::getD(spars,"DelTau",-999);
	nsample=parameter::getI(spars,"Nsample",-999);
	umax=betamax/sqrt(1.0-betamax*betamax);
	pt=fabs(p[1]);
	gammav=pt/m;
	gamma=sqrt(1.0+gammav*gammav);
	if(gammav<umax) eumin=m;
	else eumin=sqrt(1.0+umax*umax)*p[0]-umax*p[1];

	for(imc=0;imc<nsample;imc++){
		do{
			eta=etaG*randy->gauss();
			TRY_AGAIN:
			x=(1.0-2.0*randy->ran());
			y=(1.0-2.0*randy->ran());
			if(x*x+y*y>1.0) goto TRY_AGAIN;
			u[1]=umax*x;
			u[2]=umax*y;
			x=x*R;
			y=y*R;
			u[3]=sinh(eta);
			u[0]=sqrt(1.0+u[1]*u[1]+u[3]*u[3]);

			eprime=p[0]*cosh(eta);
			eu=u[0]*p[0]-u[1]*p[1];

			weight=(eprime/p[0])*exp(-(eu-eumin)/T);
			if(weight>1.0){
				printf("DISASTER! weight=%g which is > 1.0, eu=%g, eumin=%g\n",weight,eu,eumin);
				exit(1);
			}
		} while(weight<randy->ran());

		tau=GetTau(tau0,deltau);
		z=tau*sinh(eta);
		t=tau*cosh(eta);

		tt=gamma*t-gammav*x;
		x=gamma*x-gammav*t;
		t=tt;
		mclist->SetR(imc,t,x,y,z);
	}
}

double CSourceCalc_Blast::GetTau(double tau0, double deltau){
	double tau;
	do{
		tau=tau0+deltau*randy->gauss();
	}while(tau<0.0);

	return tau;
}

#endif
