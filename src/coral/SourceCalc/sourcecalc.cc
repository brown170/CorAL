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
#ifndef INCLUDE_SOURCECALC_CC
#define INCLUDE_SOURCECALC_CC

#include "constants.h"
#include "sourcecalc.h"

using namespace std;

CSourceCalc::CSourceCalc(){
	randy=new CRandom(-1234);
}

void CSourceCalc::CalcS(CCHArray *A){
	printf("i'm a dummy CalcS(CCHArray)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(C3DArray *threed){
	printf("i'm a dummy CalcS(C3DArray)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::GaussCFCalc(C3DArray *cf3d){
	printf("i'm a dummy CalcS(C3DArray *)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(int lx,int ly,int lz,CCHArray *A){
	printf("i'm a dummy CalcS(int,int,int,CCHArray *)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CalcS(CMCList *&lista,CMCList *&listb){
	printf("i'm a dummy CalcS(CMCList *a,CMCList *b)\n");
	// Dummy function to be over-written by inherited class
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A){
	A->ZeroArray();
	double rcm[4];
	int icalc=0,ncalc,icount=0,ia,ib,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
		ncalc=na*(na-1)/2;
	}
	else{
		AEQUALB=false;
		ncalc=na*nb;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<=nbmax;ib++){
			rb=listb->GetR(ib);
			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
			if(r<nrmax*delr){
				ir=int(floor(r/delr));
				if(ir<nrmax){
					A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,1.0,ir);
				}
			}
			icalc+=1;
			icount+=1;
			if(icount*10>=ncalc){
				printf("finished %4.1f percent\n",100.0*double(icalc)/double(ncalc));
				icount=0;
			}
		}
	}
	A->FillRemainderX();

	if(AEQUALB) snorm=2.0/(double(na)*double(na-1));
	else snorm=1.0/(double(na)*double(nb));

	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A,int NMC){
	A->ZeroArray();
	double rcm[4];
	int icalc=0,icount=0,ia,ib,imc,nbmax,ir,na=lista->GetNMC(),nb=listb->GetNMC(),nrmax=A->GetNRADIAL();
	double *ra,*rb;
	double r,volume,delr=A->GetRADSTEP(),snorm;
	bool AEQUALB;
	if(lista==listb){
		AEQUALB=true;
	}
	else{
		AEQUALB=false;
	}
	rcm[0]=0.0;
	printf("_______ In CombineMCLists: na=%d, nb=%d, AEQUALB=%d ______\n",na,nb,int(AEQUALB));
	A->ZeroArray();
	for(imc=0;imc<NMC;imc++){
		ia=rint(floor(randy->ran()*na));
		do{
			ib=rint(floor(randy->ran()*na));
		}while(AEQUALB && ia==ib);
		printf("ia=%d, ib=%d\n",ia,ib);
		ra=lista->GetR(ia);
		rb=listb->GetR(ib);
		rcm[1]=ra[1]-rb[1];
		rcm[2]=ra[2]-rb[2];
		rcm[3]=ra[3]-rb[3];
		r=sqrt(rcm[1]*rcm[1]+rcm[2]*rcm[2]+rcm[3]*rcm[3]);
		if(r<nrmax*delr){
			ir=int(floor(r/delr));
			if(ir<nrmax){
				A->IncrementAExpArrayFromE(rcm[1]/r,rcm[2]/r,rcm[3]/r,1.0,ir);
			}
		}
		icalc+=1;
		icount+=1;
		if(icount*10>=NMC){
			printf("finished %4.1f percent\n",100.0*double(icalc)/double(NMC));
			icount=0;

		}
	}
	A->FillRemainderX();
	snorm=1.0/double(NMC);

	for(ir=0;ir<nrmax;ir++){
		volume=(4.0*PI/3)*(pow((ir+1)*delr,3)-pow(double(ir)*delr,3));
		A->ScaleArray(snorm/volume,ir);
	}
}

void CSourceCalc::CombineMCLists(CMCList *lista,CMCList *listb,C3DArray *threed){
	threed->ZeroArray();
	double rcm[4];
	int ia,ib,nbmax,na=lista->GetNMC(),nb=listb->GetNMC();
	double *ra,*rb;
	double volume,snorm;
	bool AEQUALB=false;
	if(lista==listb) AEQUALB=true;
	rcm[0]=0.0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(AEQUALB) nbmax=ia-1;
		else nbmax=nb;
		for(ib=0;ib<nbmax;ib++){
			rb=listb->GetR(ib);
			rcm[1]=ra[1]-rb[1];
			rcm[2]=ra[2]-rb[2];
			rcm[3]=ra[3]-rb[3];
			threed->IncrementElement(rcm[1],rcm[2],rcm[3],1.0);	
		}
		if(10*(ia+1)%(10*int(na/10))==0)
			printf("finished %g percent\n",100*double(ia+1)/(10*int(na/10)));
	}

	if(AEQUALB) snorm=2.0/double(na*(na-1));
	else snorm=1.0/double(na*nb);

	volume=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
	threed->ScaleArray(snorm/volume);
}

void CSourceCalc::ReadSPars(char *sparsfilename){
	parameter::ReadParsFromFile(spars,sparsfilename);
}

void CSourceCalc::NormCheck(CCHArray *A){
	double check,delr,volume;
	int ir,NRMAX;
	NRMAX=A->GetNRADIAL();
	delr=A->GetRADSTEP();
	check=0;
	for(ir=0;ir<NRMAX;ir++){
		volume=4.0*PI*(pow(delr*(ir+1.0),3)-pow(delr*ir,3))/3.0;
		check+=volume*A->GetElement(0,0,0,ir);
	}
	printf("normalization check = %g\n",check);
}

void CSourceCalc::CalcEffGaussPars(CCHArray *A){
	double Rx,Ry,Rz,Xoff,Yoff,Zoff;
	CalcEffGaussPars(A,Rx,Ry,Rz,Xoff,Yoff,Zoff);
}

void CSourceCalc::CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,double &Rz,double &Xoff,double &Yoff,double &Zoff){
	double xbar,ybar,zbar,x2bar,y2bar,z2bar,r2bar,r,r2,r3,r4,DELR;
	int ir,NRMAX;
	bool XSYM,YSYM,ZSYM;
	NRMAX=A->GetNRADIAL();
	DELR=A->GetRADSTEP();
	XSYM=A->GetXSYM();
	YSYM=A->GetYSYM();
	ZSYM=A->GetZSYM();
	const double PI=4.0*atan(1.0);
	xbar=ybar=zbar=x2bar=y2bar=z2bar=r2bar=0.0;
	for(ir=0;ir<NRMAX;ir++){
		r=(0.5+ir)*DELR;
		r2=r*r; 
		r3=r2*r;
		r4=r2*r2;
		if(!XSYM) xbar+=r3*A->GetElement(1,0,0,ir);
		if(!YSYM) ybar+=r3*A->GetElement(0,1,0,ir);
		if(!ZSYM) zbar+=r3*A->GetElement(0,0,1,ir);
		x2bar+=r4*A->GetElement(2,0,0,ir);
		y2bar+=r4*A->GetElement(0,2,0,ir);
		z2bar+=r4*A->GetElement(0,0,2,ir);
		r2bar+=r4*A->GetElement(0,0,0,ir);
	}
	xbar*=4.0*PI*DELR/3.0;
	ybar*=4.0*PI*DELR/3.0;
	zbar*=4.0*PI*DELR/3.0;
	printf("__________  EFFECTIVE GAUSSIAN PARAMETERS ____________\n");
	printf("Rinv=%g\n",sqrt(2.0*PI*DELR*r2bar/3.0));
	x2bar=4.0*PI*DELR*(2.0*x2bar/15.0+(r2bar/3.0))-xbar*xbar;
	y2bar=4.0*PI*DELR*(2.0*y2bar/15.0+(r2bar/3.0))-ybar*ybar;
	z2bar=4.0*PI*DELR*(2.0*z2bar/15.0+(r2bar/3.0))-zbar*zbar;
	printf("Gaussian distribution with same offsets and 1-part. radii\n");
	printf("offset_xyz=(%g,%g,%g), R_xyz=(%g,%g,%g)\n",
		xbar,ybar,zbar,
		sqrt(fabs(0.5*x2bar)),sqrt(fabs(0.5*y2bar)),sqrt(fabs(0.5*z2bar)));
	printf("______________________________________________________\n");

	Xoff=xbar;
	Yoff=ybar;
	Zoff=zbar;
	Rx=sqrt(fabs(0.5*x2bar));
	Ry=sqrt(fabs(0.5*y2bar));
	Rz=sqrt(fabs(0.5*z2bar));

}

void CSourceCalc::NormCheck(C3DArray *threed){
	int nsx,nsy,nsz,isx,isy,isz,ix,iy,iz;
	int nxmax=threed->GetNXMAX();
	int nymax=threed->GetNYMAX();
	int nzmax=threed->GetNZMAX();
	double prefactor=threed->GetDELX()*threed->GetDELY()*threed->GetDELZ();
	double norm=0.0;
	nsx=nsy=nsz=2;
	if(threed->GetXSYM()) nsx=1;
	if(threed->GetYSYM()) nsy=1;
	if(threed->GetZSYM()) nsz=1;
	if(nsx==1) prefactor*=2.0;
	if(nsy==1) prefactor*=2.0;
	if(nsz==1) prefactor*=2.0;
	for(ix=0;ix<nxmax;ix++){
		for(iy=0;iy<nymax;iy++){
			for(iz=0;iz<nzmax;iz++){
				for(isz=0;isz<nsz;isz++){
					for(isy=0;isy<nsy;isy++){
						for(isx=0;isx<nsx;isx++){
							norm+=threed->GetElement(isx,ix,isy,iy,isz,iz)*prefactor;
						}
					}
				}
			}
		}
	}
	printf("Norm Check of 3DArray = %g\n",norm);
}

#endif
