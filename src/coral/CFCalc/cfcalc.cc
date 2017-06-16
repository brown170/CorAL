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
#ifndef __INCLUDE_CFCalc_CC__
#define __INCLUDE_CFCalc_CC__

#include "cfcalc.h"
using namespace std;

double CFCalc::GetChiSquared(C3DArray *CFexp,C3DArray *Error, C3DArray *CFtheory){
	if(!ArrayCalc::CompareArrayParameters(CFexp,Error)) exit(1);
	if(!ArrayCalc::CompareArrayParameters(CFexp,CFtheory)) exit(1);
	double chisquared=0.0,numer,denom;
	int ix,iy,iz,isx,isy,isz,nsx,nsy,nsz,nxmax,nymax,nzmax;
	int ndof=0;
	nsx=nsy=nsz=2;
	if(CFexp->GetXSYM()) nsx=1;
	if(CFexp->GetYSYM()) nsy=1;
	if(CFexp->GetZSYM()) nsz=1;
	nxmax=CFexp->GetNXMAX();
	nymax=CFexp->GetNYMAX();
	nzmax=CFexp->GetNZMAX();
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							numer=CFexp->GetElement(isx,ix,isy,iy,isz,iz)
								-CFtheory->GetElement(isx,ix,isy,iy,isz,iz);
							numer=numer*numer;
							denom=Error->GetElement(isx,ix,isy,iy,isz,iz);
							if(denom<-1.0E-8) printf("CFCalc::GetChiSquared : \nSuspiciously small or negative error = %g, isx=%d,ix=%d, isy=%d,iy=%d, isz=%d,iz=%d\n",
								denom,isx,ix,isy,iy,isz,iz);
							denom=denom*denom;
							if(numer==numer && denom==denom){
								chisquared+=numer/denom;
							 	ndof+=1;
							}
						}
					}
				}
			}
		}
	}
	printf("chi^2/Ndof=%g\n",chisquared/double(ndof));
	return chisquared;
}

double CFCalc::GetChiSquared(int lx,int ly,int lz,CCHArray *CFexp, CCHArray *Error,CCHArray *CFtheory){
	int ir,nrmax;
	double chisquared=0.0,numer,denom;
	nrmax=CFexp->GetNRADIAL();
	for(ir=0;ir<nrmax;ir++){
		denom=Error->GetElement(lx,ly,lz,ir);
		if(denom<1.0E-8) printf("CFCalc::GetChiSquared : \nSuspiciously small or negative error = %g, l=(%d,%d,%d)d\n",denom,lx,ly,lz);
		numer=CFexp->GetElement(lx,ly,lz,ir)
			-CFtheory->GetElement(lx,ly,lz,ir);
		chisquared+=numer*numer/(denom*denom);
	}
	printf("chi^2/Ndof=%g\n",chisquared/double(nrmax));
	return chisquared;
}

double CFCalc::GetChiSquared(CCHArray *cexp,CCHArray *error,CCHArray *ctheory){
	int L,Lmax,lx,ly,lz,lxprime,lyprime,lzprime,lxprime0,lyprime0,iq,nqmax;
	double lfact,lprimefact;
	int dlx,dly,dlz;
	dlx=dly=dlz=1;
	if(cexp->GetXSYM()) dlx=2;
	if(cexp->GetYSYM()) dly=2;
	if(cexp->GetZSYM()) dlz=2;
	double overlap,Fvalue,Gvalue,denom,prefact,q;
	nqmax=error->GetNRADIAL();
	Lmax=ctheory->GetLMAX();
	CCHArray *F=new CCHArray(Lmax,nqmax,cexp->GetRADSTEP(),cexp->GetXSYM(),cexp->GetYSYM(),cexp->GetZSYM());
	ArrayCalc::SubtractArrays(ctheory,cexp,F);
	overlap=0.0;
	for(iq=0;iq<nqmax;iq++){
		denom=error->GetElement(0,0,0,iq);
		q=F->GetRADSTEP()*(iq+0.5);
		prefact=4.0*PI*q*q*F->GetRADSTEP()/denom;
		for(L=0;L<=Lmax;L++){
			for(lx=0;lx<=L;lx+=dlx){
				for(ly=0;ly<=L-lx;ly+=dly){
					lz=L-lx-ly;
					if(lz%dlz==0){
						lfact=F->chcalc->Factorial(L)
							/(F->chcalc->Factorial(lx)*F->chcalc->Factorial(ly)*F->chcalc->Factorial(lz));
						Fvalue=F->GetElement(lx,ly,lz,iq);
						lxprime0=0;
						if(lx%2==1) lxprime0=1;
						lyprime0=0;
						if(ly%2==1) lyprime0=1;
						for(lxprime=lxprime0;lxprime<=L;lxprime+=2){
							for(lyprime=lyprime0;lyprime<=L-lxprime;lyprime+=2){
								lzprime=L-lxprime-lyprime;
								lprimefact=F->chcalc->Factorial(L)
									/(F->chcalc->Factorial(lxprime)*F->chcalc->Factorial(lyprime)*F->chcalc->Factorial(lzprime));
								Gvalue=F->GetElement(lxprime,lyprime,lzprime,iq);
								overlap+=prefact*lfact*lprimefact*Fvalue*Gvalue*F->chcalc->GetOverlap(lx,ly,lz,lxprime,lyprime,lzprime);
							}
						}
					}
				}
			}
		}
	}
	printf("chi^2=%g\n",overlap);
	delete(F);
	return overlap;
}

#endif
