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
#ifndef __INCLUDE_S2C_CC
#define __INCLUDE_S2C_CC

#include "source2cf.h"

using namespace std;

void S2CF::s2c(int lx,int ly,int lz,CCHArray *S,CKernel *kernel,CCHArray *CF){
	int lmax=kernel->GetLMAX();
	if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
		printf("FATAL: Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
			lmax,S->GetLMAX(),CF->GetLMAX());
		exit(1);
	}
	if( (CF->GetLMAX()!=lmax) ){
		printf("WARNING: Array parameters for kernel calculations don't match\n");
		printf("___ CORRELATION FUNCTION PARAMETERS ___\n");
		CF->PrintPars();
		printf("_____ SOURCE FUNCTION PARAMETERS _____\n");
		S->PrintPars();
		printf("For kernel, LMAX=%d\n",kernel->GetLMAX());
	}
	if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
		printf("FATAL: kernel has no odd L components, but CF wants them\n");
		printf("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");
		exit(1);
	}
	int iq,ir,nqmax,nrmax,L;
	double r,delr,q,delq,norm;
	const double PI=4.0*atan(1.0);
	bool match=0;
	delr=S->GetRADSTEP();
	nrmax=S->GetNRADIAL();
	delq=CF->GetRADSTEP();
	nqmax=CF->GetNRADIAL();
	if(fabs(delr-kernel->GetDELR())<1.0E-5 
		&& fabs(delq-kernel->GetDELQ())<1.0E-5
		&& nrmax==kernel->GetNRMAX()
		&& nqmax==kernel->GetNQMAX()) match=1;
	CF->ZeroArray(lx,ly,lz);
	L=lx+ly+lz;
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(ir=0;ir<nrmax;ir++){
			r=(0.5+ir)*delr;
			norm=4.0*PI*r*r*delr;
			if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)*norm*S->GetElement(lx,ly,lz,ir));
			else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)*norm*S->GetElement(lx,ly,lz,r));
		}
	}
}

void S2CF::s2c(CCHArray *S,CKernel *kernel,CCHArray *CF){
	int lmax=kernel->GetLMAX();
	if(lmax>CF->GetLMAX() || lmax>S->GetLMAX()){
		printf("FATAL: Kernel LMAX=%d larger than either S LMAX=%d or CF LMAX=%d\n",
			lmax,S->GetLMAX(),CF->GetLMAX());
		exit(1);
	}
	if( (CF->GetLMAX()!=lmax) ){
		printf("WARNING: Array parameters for kernel calculations don't match\n");
		printf("___ CORRELATION FUNCTION PARAMETERS ___\n");
		CF->PrintPars();
		printf("_____ SOURCE FUNCTION PARAMETERS _____\n");
		S->PrintPars();
		printf("For kernel, LMAX=%d\n",kernel->GetLMAX());
	}
	if(kernel->GetIDENTICAL() && (!CF->GetXSYM() || !CF->GetYSYM() || !CF->GetZSYM())){
		printf("FATAL: kernel has no odd L components, but CF wants them\n");
		printf("Make sure CF array has XSYM=YSYM=ZSYM='true'\n");
		exit(1);
	}
	int iq,ir,nqmax,nrmax,L,lx,ly,lz,dlx,dly,dlz;
	double r,delr,q,delq,norm;
	const double PI=4.0*atan(1.0);
	bool match=0;
	delr=S->GetRADSTEP();
	nrmax=S->GetNRADIAL();
	delq=CF->GetRADSTEP();
	nqmax=CF->GetNRADIAL();
	dlx=dly=dlz=1;
	if(CF->GetXSYM()) dlx=2;
	if(CF->GetYSYM()) dly=2;
	if(CF->GetZSYM()) dlz=2;
	if(fabs(delr-kernel->GetDELR())<1.0E-5 
		&& fabs(delq-kernel->GetDELQ())<1.0E-5
		&& nrmax==kernel->GetNRMAX()
		&& nqmax==kernel->GetNQMAX()) match=1;
	CF->ZeroArray();
	for(iq=0;iq<nqmax;iq++){
		q=(0.5+iq)*delq;
		for(lx=0;lx<2;lx+=dlx){
			for(ly=0;ly<=lmax-lx;ly+=dly){
				for(lz=0;lz<=lmax-lx-ly;lz+=dlz){
					L=lx+ly+lz;
					for(ir=0;ir<nrmax;ir++){
						r=(0.5+ir)*delr;
						norm=4.0*PI*r*r*delr;
						if(match) CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,iq,ir)
							*norm*S->GetElement(lx,ly,lz,ir));
						else CF->IncrementElement(lx,ly,lz,iq,kernel->GetValue(L,q,r)
							*norm*S->GetElement(lx,ly,lz,r));
					}
				}
			}
		}
	}
}

void S2CF::s2c(C3DArray *S,CKernelWF *kernel,C3DArray *CF){
	int ix,iy,iz,isx,isy,isz,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	long long int icalc,ncalc;
	int nxmax=S->GetNXMAX();
	int nymax=S->GetNYMAX();
	int nzmax=S->GetNZMAX();
	double delx=S->GetDELX();
	double dely=S->GetDELY();
	double delz=S->GetDELZ();
	int nqxmax=CF->GetNXMAX();
	int nqymax=CF->GetNYMAX();
	int nqzmax=CF->GetNZMAX();
	double delqx=CF->GetDELX();
	double delqy=CF->GetDELY();
	double delqz=CF->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,svalue;
	bool IDENTICAL=kernel->GetIDENTICAL();
	bool XSYM=S->GetXSYM();
	bool YSYM=S->GetYSYM();
	bool ZSYM=S->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		printf("S2CF::s2c, kernel says particles are identical, but symmetries are not all even\n");
		printf("XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		exit(1);
	}

	norm=delx*dely*delz;
	nsx=nsy=nsz=2;
	if(XSYM){
		nsx=1;
		norm*=2.0;
	}
	if(YSYM){
		nsy=1;
		norm*=2.0;
	}
	if(ZSYM){
		nsz=1;
		norm*=2.0;
	}
	CF->ZeroArray();
	ncalc=nxmax*nymax*nzmax*nsx*nsy*nsz;
	ncalc=ncalc/10;
	icalc=0;

	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			x=(0.5+ix)*delx;
			if(isx==1) x=-x;

			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					y=(0.5+iy)*dely;
					if(isy==1) y=-y;

					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							z=(0.5+iz)*delz;
							if(isz==1) z=-z;

							r=sqrt(x*x+y*y+z*z);
							svalue=S->GetElement(isx,ix,isy,iy,isz,iz);
							//svalue=exp(-x*x/(4.0*9.0)-y*y/(4.0*25)-z*z/(4.0*49));
							//svalue=svalue/(pow(4.0*PI*9.0*25.0*49.0);

							for(jsx=0;jsx<nsx;jsx++){
								for(jx=0;jx<nqxmax;jx++){
									qx=(0.5+jx)*delqx;
									if(jsx==1) qx=-qx;

									for(jsy=0;jsy<nsy;jsy++){
										for(jy=0;jy<nqymax;jy++){
											qy=(0.5+jy)*delqy;
											if(jsy==1) qy=-qy;

											for(jsz=0;jsz<nsz;jsz++){
												for(jz=0;jz<nqzmax;jz++){
													qz=(0.5+jz)*delqz;
													if(jsz==1) qz=-qz;
													q=sqrt(qx*qx+qy*qy+qz*qz);
													wf2=0.0;

													if(XSYM&YSYM&&ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);																																																				
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
														if(!IDENTICAL) wf2=0.5*wf2;
													}
													else if(!XSYM && YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
													}
													else{
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=kernel->GetPsiSquared(q,r,ctheta);
													}
													CF->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
														norm*(wf2-1.0)*svalue);
												}
											}
										}
									}
								}
							}
							icalc+=1;
							if(icalc%ncalc==0) printf("s2c, finished %g percent\n",10*double(icalc)/double(ncalc));
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c(C3DArray *S,CWaveFunction *wf,C3DArray *CF){
	int ix,iy,iz,isx,isy,isz,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nxmax=S->GetNXMAX();
	int nymax=S->GetNYMAX();
	int nzmax=S->GetNZMAX();
	double delx=S->GetDELX();
	double dely=S->GetDELY();
	double delz=S->GetDELZ();
	int nqxmax=CF->GetNXMAX();
	int nqymax=CF->GetNYMAX();
	int nqzmax=CF->GetNZMAX();
	double delqx=CF->GetDELX();
	double delqy=CF->GetDELY();
	double delqz=CF->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,svalue;
	bool IDENTICAL=wf->GetIDENTICAL();
	bool XSYM=S->GetXSYM();
	bool YSYM=S->GetYSYM();
	bool ZSYM=S->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		printf("S2CF::s2c, wf says particles are identical, but symmetries are not all even\n");
		printf("XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		exit(1);
	}

	norm=delx*dely*delz;
	nsx=nsy=nsz=2;
	if(XSYM){
		nsx=1;
		norm*=2.0;
	}
	if(YSYM){
		nsy=1;
		norm*=2.0;
	}
	if(ZSYM){
		nsz=1;
		norm*=2.0;
	}
	CF->ZeroArray();

	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			x=(0.5+ix)*delx;
			if(isx==1) x=-x;

			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					y=(0.5+iy)*dely;
					if(isy==1) y=-y;

					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							z=(0.5+iz)*delz;
							if(isz==1) z=-z;

							r=sqrt(x*x+y*y+z*z);
							//svalue=S->GetElement(isx,ix,isy,iy,isz,iz);
							svalue=exp(-0.25*((x*x/9.0)+(y*y/25.0)+(z*z/49.0)));
							svalue=svalue/(3.0*5.0*7.0*pow(4.0*PI,1.5));

							for(jsx=0;jsx<nsx;jsx++){
								for(jx=0;jx<nqxmax;jx++){
									qx=(0.5+jx)*delqx;
									if(jsx==1) qx=-qx;

									for(jsy=0;jsy<nsy;jsy++){
										for(jy=0;jy<nqymax;jy++){
											qy=(0.5+jy)*delqy;
											if(jsy==1) qy=-qy;

											for(jsz=0;jsz<nsz;jsz++){
												for(jz=0;jz<nqzmax;jz++){
													qz=(0.5+jz)*delqz;
													if(jsz==1) qz=-qz;
													q=sqrt(qx*qx+qy*qy+qz*qz);
													wf2=0.0;

													if(XSYM&YSYM&&ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);																																																				
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														if(!IDENTICAL) wf2+=0.25*wf->GetPsiSquared(q,r,-ctheta);
														if(!IDENTICAL) wf2=0.5*wf2;
													}
													else if(!XSYM && YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.25*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(XSYM && !YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(-qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && YSYM && !ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x-qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else if(!XSYM && !YSYM && ZSYM){
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
														ctheta=(qx*x+qy*y-qz*z)/(q*r);
														wf2+=0.5*wf->GetPsiSquared(q,r,ctheta);
													}
													else{
														ctheta=(qx*x+qy*y+qz*z)/(q*r);
														wf2+=wf->GetPsiSquared(q,r,ctheta);
													}
													CF->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
														norm*(wf2-1.0)*svalue);
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf){
	bool samelists;
	int ia,ib,na,nb,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	int npairs,ipair,jpair;
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetNRMAX()*kernel->GetDELR();
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !kernel->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	na=lista->GetNMC();
	nb=lista->GetNMC();
	ipair=jpair=0;
	//printf("begin loop in s2c(CMCList *,CMCList *,CKernelWF *,C3DArray *)\n");
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			if(r<rmax){
				for(jsx=0;jsx<nsx;jsx++){
					for(jx=0;jx<nqxmax;jx++){
						qx=(0.5+jx)*delqx;
						if(jsx==1) qx=-qx;
						for(jsy=0;jsy<nsy;jsy++){
							for(jy=0;jy<nqymax;jy++){
								qy=(0.5+jy)*delqy;
								if(jsy==1) qy=-qy;
								for(jsz=0;jsz<nsz;jsz++){
									for(jz=0;jz<nqzmax;jz++){
										qz=(0.5+jz)*delqz;
										if(jsz==1) qz=-qz;
										q=sqrt(qx*qx+qy*qy+qz*qz);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										wf2=kernel->GetPsiSquared(q,r,ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										if(samelists && !kernel->GetIDENTICAL()){
											wf2=kernel->GetPsiSquared(q,r,-ctheta);
											cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										}

									}
								}
							}
						}
					}
				}
				ipair+=1;
				if(ipair>=npairs/10){
					jpair+=1;
					printf("Finished %d percent\n",jpair*10);
					ipair=0;
				}
			}
		}
	}
	//printf("ending big loop in s2c\n");
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernel *kernel,C3DArray *cf){
	bool samelists;
	int ia,ib,na,nb,npairs,jx,jy,jz,jsx,jsy,jsz,ipair,jpair;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetDELR()*kernel->GetNRMAX();
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !kernel->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	//printf("beginning big loop in s2c\n");
	na=lista->GetNMC();
	nb=listb->GetNMC();
	ipair=jpair=0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			if(r<rmax){
				for(jsx=0;jsx<nsx;jsx++){
					for(jx=0;jx<nqxmax;jx++){
						qx=(0.5+jx)*delqx;
						if(jsx==1) qx=-qx;
						for(jsy=0;jsy<nsy;jsy++){
							for(jy=0;jy<nqymax;jy++){
								qy=(0.5+jy)*delqy;
								if(jsy==1) qy=-qy;
								for(jsz=0;jsz<nsz;jsz++){
									for(jz=0;jz<nqzmax;jz++){
										qz=(0.5+jz)*delqz;
										if(jsz==1) qz=-qz;
										q=sqrt(qx*qx+qy*qy+qz*qz);
										ctheta=(qx*x+qy*y+qz*z)/(q*r);
										wf2=kernel->GetPsiSquared(q,r,ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										if(samelists && !kernel->GetIDENTICAL()){
											wf2=kernel->GetPsiSquared(q,r,-ctheta);
											cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
										}

									}
								}
							}
						}
					}
				}
			}
			ipair+=1;
			if(ipair==npairs/10){
				jpair+=1;
				printf("Finished %d percent\n",jpair*10);
				ipair=0;
			}
		}
	}
	//printf("ending big loop in s2c\n");
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf,int NMC){
	CRandom randy(-time(NULL));
	bool samelists;
	int ia,ib,na,nb,jx,jy,jz,jsx,jsy,jsz,imc,jmc,npairs;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	double rmax=kernel->GetDELR()*kernel->GetNRMAX();
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(NMC);

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();
	na=lista->GetNMC();
	nb=listb->GetNMC();

	for(jsx=0;jsx<nsx;jsx++){
		for(jx=0;jx<nqxmax;jx++){
			qx=(0.5+jx)*delqx;
			if(jsx==1) qx=-qx;
			for(jsy=0;jsy<nsy;jsy++){
				for(jy=0;jy<nqymax;jy++){
					qy=(0.5+jy)*delqy;
					if(jsy==1) qy=-qy;
					for(jsz=0;jsz<nsz;jsz++){
						for(jz=0;jz<nqzmax;jz++){
							qz=(0.5+jz)*delqz;
							if(jsz==1) qz=-qz;
							q=sqrt(qx*qx+qy*qy+qz*qz);
							imc=jmc=0;
							for(imc=0;imc<NMC;imc++){
								jmc+=1;
								ia=rint(floor(randy.ran()*na));
								do{
									ib=rint(floor(randy.ran()*na));
								} while(ia==ib && samelists);
								ra=lista->GetR(ia);
								rb=listb->GetR(ib);
								x=ra[1]-rb[1];
								y=ra[2]-rb[2];
								z=ra[3]-rb[3];
								r=sqrt(x*x+y*y+z*z);
								if(r<rmax){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2=kernel->GetPsiSquared(q,r,ctheta);
									cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									if(samelists && !kernel->GetIDENTICAL()){
										wf2=kernel->GetPsiSquared(q,r,-ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									}

								}
							}
						}
					}
				}
			}
		}
		if(jmc==NMC/10){
			printf("finished %g percent\n",double(imc+1)*100/double(NMC));
			jmc=0;
		}
	}
}

void S2CF::s2c(CMCList *lista,CMCList *listb,CWaveFunction *wf,C3DArray *cf){
	bool samelists;
	int ia,ib,na,nb,npairs,ipair,jpair,jx,jy,jz,jsx,jsy,jsz;
	int nsx,nsy,nsz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double x,y,z,qx,qy,qz,norm,q,r,ctheta,wf2,*ra,*rb;
	if(lista==listb){
		npairs=lista->GetNMC()*(lista->GetNMC()-1)/2;
		samelists=true;
	}
	else{
		npairs=lista->GetNMC()*listb->GetNMC();
		samelists=false;
	}
	norm=lista->GetNorm()*listb->GetNorm()/double(npairs);
	if(samelists && !wf->GetIDENTICAL()) norm=0.5*norm;

	nsx=nsy=nsz=2;
	if(cf->GetXSYM()){
		nsx=1;
	}
	if(cf->GetYSYM()){
		nsy=1;
	}
	if(cf->GetZSYM()){
		nsz=1;
	}
	cf->ZeroArray();

	//printf("begin loop in s2c(CMCList *,CMCList *,CWaveFunction *,C3DArray *)\n");
	na=lista->GetNMC();
	nb=listb->GetNMC();
	ipair=jpair=0;
	for(ia=0;ia<na;ia++){
		ra=lista->GetR(ia);
		if(samelists) nb=ia-1;
		for(ib=0;ib<nb;ib++){
			rb=listb->GetR(ib);
			x=ra[1]-rb[1];
			y=ra[2]-rb[2];
			z=ra[3]-rb[3];
			r=sqrt(x*x+y*y+z*z);
			for(jsx=0;jsx<nsx;jsx++){
				for(jx=0;jx<nqxmax;jx++){
					qx=(0.5+jx)*delqx;
					if(jsx==1) qx=-qx;
					for(jsy=0;jsy<nsy;jsy++){
						for(jy=0;jy<nqymax;jy++){
							qy=(0.5+jy)*delqy;
							if(jsy==1) qy=-qy;
							for(jsz=0;jsz<nsz;jsz++){
								for(jz=0;jz<nqzmax;jz++){
									qz=(0.5+jz)*delqz;
									if(jsz==1) qz=-qz;
									q=sqrt(qx*qx+qy*qy+qz*qz);
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2=wf->GetPsiSquared(q,r,ctheta);
									cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									if(samelists && !wf->GetIDENTICAL()){
										wf2=wf->GetPsiSquared(q,r,-ctheta);
										cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,norm*(wf2-1.0));
									}

								}
							}
						}
					}
				}
			}
			ipair+=1;
			if(ipair>=npairs/10){
				jpair+=1;
				printf("Finished %d percent\n",jpair*10);
				ipair=0;
			}
		}
	}
	//for(jx=0;jx<nqxmax;jx++) printf("C(iqx=%d)=%g\n",jx,cf->GetElement(0,jx,0,0,0,0));
	//printf("ending big loop in s2c\n");
}

void S2CF::s2c_bosons(CMCList *list,C3DArray *cf){
	int i,n,jx,jy,jz;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double c,r2,norm,arg,qx,qy,qz,*r;
	complex<double> alpha,ci(0.0,1.0);
	n=list->GetNMC();
	norm=list->GetNorm()*list->GetNorm();

	cf->ZeroArray();
	for(jx=0;jx<nqxmax;jx++){
		qx=(0.5+jx)*delqx;
		for(jy=0;jy<nqymax;jy++){
			qy=(0.5+jy)*delqy;
			for(jz=0;jz<nqzmax;jz++){
				qz=(0.5+jz)*delqz;
				alpha=0.0;
				for(i=0;i<n;i++){
					r=list->GetR(i);
					r2=r[1]*r[1]+r[2]*r[2]+r[3]*r[3];
					arg=2.0*(r[1]*qx+r[2]*qy+r[3]*qz)/HBARC;
					if(arg>2.0*PI) arg=arg-2.0*PI*floor(arg/(2.0*PI));
					if(arg<0.0) arg=arg+2.0*PI*floor(fabs(arg/(2.0*PI)));
					if(r2<1.0E10)	alpha+=exp(ci*arg);
				}
				c=abs(alpha*conj(alpha))-double(n);
				c=c/(double(n)*double(n-1));
				cf->IncrementElement(0,jx,0,jy,0,jz,norm*c);
			}
		}
	}
}

void S2CF::s2c_gauss(CSourceCalc *sourcecalc,CKernelWF *kernel,C3DArray *cf){
	int imc;
	double root2=sqrt(2.0);
	double Rx=parameter::getD(sourcecalc->spars,"Rx",4.0);
	double Ry=parameter::getD(sourcecalc->spars,"Ry",4.0);
	double Rz=parameter::getD(sourcecalc->spars,"Rz",4.0);
	int NMC=parameter::getI(sourcecalc->spars,"NMC",500);
	double lambda=parameter::getD(sourcecalc->spars,"lambda",1.0);
	double Xoff=parameter::getD(sourcecalc->spars,"Xoff",0.0);
	double Yoff=parameter::getD(sourcecalc->spars,"Yoff",0.0);
	double Zoff=parameter::getD(sourcecalc->spars,"Zoff",0.0);
	//double Euler_Phi=parameter::set(spars,"Euler_Phi",0.0);
	//double Euler_Theta=parameter::set(spars,"Euler_Theta",0.0);
	//double Euler_Psi=parameter::set(spars,"Euler_Psi",0.0);
	CRandom *randy=sourcecalc->randy;
	double x,y,z,r,q,ctheta,qx,qy,qz,wf2,xarray[3],gauss[2];
	int igauss=2,ix;

	bool IDENTICAL=kernel->GetIDENTICAL();
	bool XSYM=cf->GetXSYM();
	bool YSYM=cf->GetYSYM();
	bool ZSYM=cf->GetZSYM();
	if(IDENTICAL&&!(XSYM&&YSYM&&ZSYM)){
		printf("S2CF::s2c, kernel says particles are identical, but symmetries are not all even\n");
		printf("XSYM=%d, YSYM=%d, ZSYM=%d\n",int(XSYM),int(YSYM),int(ZSYM));
		exit(1);
	}
	int jsx,jsy,jsz,jx,jy,jz;
	int nsx=2,nsy=2,nsz=2;
	if(XSYM) nsx=1;
	if(YSYM) nsy=1;
	if(ZSYM) nsz=1;
	int nqxmax=cf->GetNXMAX();
	int nqymax=cf->GetNYMAX();
	int nqzmax=cf->GetNZMAX();
	double delqx=cf->GetDELX();
	double delqy=cf->GetDELY();
	double delqz=cf->GetDELZ();
	double norm=lambda/double(NMC);
	cf->ZeroArray();

	for(jsx=0;jsx<nsx;jsx++){
		for(jx=0;jx<nqxmax;jx++){
			qx=(0.5+jx)*delqx;
			if(jsx==1) qx=-qx;

			for(jsy=0;jsy<nsy;jsy++){
				for(jy=0;jy<nqymax;jy++){
					qy=(0.5+jy)*delqy;
					if(jsy==1) qy=-qy;

					for(jsz=0;jsz<nsz;jsz++){
						for(jz=0;jz<nqzmax;jz++){
							qz=(0.5+jz)*delqz;
							if(jsz==1) qz=-qz;
							q=sqrt(qx*qx+qy*qy+qz*qz);

							for(imc=0;imc<NMC;imc++){
								for(ix=0;ix<3;ix++){
									if(igauss==2){
										randy->gauss2(&gauss[0],&gauss[1]);
										igauss=0;
									}
									xarray[ix]=gauss[igauss]; igauss++;
								}
								x=root2*Rx*xarray[0]+Xoff;
								y=root2*Ry*xarray[1]+Yoff;
								z=root2*Rz*xarray[2]+Zoff;
								r=sqrt(x*x+y*y+z*z);
								wf2=0.0;

								if(XSYM&YSYM&&ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);																																																				
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									if(!IDENTICAL) wf2+=0.25*kernel->GetPsiSquared(q,r,-ctheta);
									if(!IDENTICAL) wf2=0.5*wf2;
								}
								else if(!XSYM && YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && !YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.25*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(XSYM && !YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(-qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(!XSYM && YSYM && !ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x-qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else if(!XSYM && !YSYM && ZSYM){
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
									ctheta=(qx*x+qy*y-qz*z)/(q*r);
									wf2+=0.5*kernel->GetPsiSquared(q,r,ctheta);
								}
								else{
									ctheta=(qx*x+qy*y+qz*z)/(q*r);
									wf2+=kernel->GetPsiSquared(q,r,ctheta);
								}
								cf->IncrementElement(jsx,jx,jsy,jy,jsz,jz,
									norm*(wf2-1.0));
							}
						}
					}
				}
			}
		}
	}
}

void S2CF::s2c_bowlersinyukov(CSourceCalc *sourcecalc,CKernel *kernel,C3DArray *cf3d){
	double Rx=parameter::getD(sourcecalc->spars,"Rx",4.0);
	double Ry=parameter::getD(sourcecalc->spars,"Ry",4.0);
	double Rz=parameter::getD(sourcecalc->spars,"Rz",4.0);
	double HBARC2=HBARC*HBARC;
	double lambda=parameter::getD(sourcecalc->spars,"lambda",1.0);
	//double Rbar=pow(Rx*Ry*Rz,1.0/3.0);
	double Rbar=5.0;
	double snorm=1.0/pow(4.0*PI*Rbar*Rbar,1.5);
	double nrmax=kernel->GetNRMAX();
	double nqmax=kernel->GetNQMAX();
	double delr=kernel->GetDELR();
	double delq=kernel->GetDELQ();
	int nqxmax=cf3d->GetNXMAX();
	int nqymax=cf3d->GetNYMAX();
	int nqzmax=cf3d->GetNZMAX();
	double delqx=cf3d->GetDELX();
	double delqy=cf3d->GetDELY();
	double delqz=cf3d->GetDELZ();
	double r,ss;
	CCHArray *sf=new CCHArray(0,nrmax,delr,true,true,true);
	double norm=0.0;
	for(int ir=0;ir<nrmax;ir++){
		r=(0.5+ir)*delr;
		ss=snorm*exp(-0.25*r*r/(Rbar*Rbar));
		sf->SetElement(0,0,0,ir,ss);
		norm+=ss*4.0*PI*r*r*delr;
	}
	CCHArray *cf=new CCHArray(0,nqmax,delq,true,true,true);
	S2CF::s2c(sf,kernel,cf);
	//cf->Print(0,0,0);
	//Misc::Pause();

	int jx,jy,jz,iq,iq1,iq2;
	double qx,qy,qz,q,cfbar,w1,w2,q1,q2,bowlsin;	
	for(jx=0;jx<nqxmax;jx++){
		qx=(0.5+jx)*delqx;
		for(jy=0;jy<nqymax;jy++){
			qy=(0.5+jy)*delqy;
			for(jz=0;jz<nqzmax;jz++){
				qz=(0.5+jz)*delqz;
				q=sqrt(qx*qx+qy*qy+qz*qz);
				iq=lrint(q/delq);
				iq1=iq-1;
				iq2=iq;
				if(iq1<0){
					iq1+=1;
					iq2+=1;
				}
				q1=(0.5+iq1)*delq;
				q2=(0.5+iq2)*delq;
				w1=(q2-q)/delq;
				w2=1.0-w1;
				cfbar=w1*cf->GetElement(0,0,0,iq1)+w2*cf->GetElement(0,0,0,iq2);
				bowlsin=lambda*(1.0+exp(-4.0*qx*qx*Rx*Rx/HBARC2-4.0*qy*qy*Ry*Ry/HBARC2-4.0* qz*qz*Rz*Rz/HBARC2))*(1.0+cfbar);
				bowlsin-=lambda;
				cf3d->SetElement(0,jx,0,jy,0,jz,bowlsin);
			}
		}
	}
	//printf("lambda=%g\n",lambda);
	//cf3d->PrintProjections();
	//Misc::Pause();
	delete(sf);
	delete(cf);
}

#endif
