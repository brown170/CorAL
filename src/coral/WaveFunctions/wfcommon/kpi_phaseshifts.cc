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
#ifndef __CWAVEFUNCTION_WF_KPI_PHASESHIFT_CC__
#define __CWAVEFUNCTION_WF_KPI_PHASESHIFT_CC__
#include "wavefunction.h"

void WaveFunctionRoutines::getphaseshift_kpi(int twoI,int ell,double q,double *delta,double *ddeltadq){
  double a,reff;
	const double R=(1000.0/4.3),GammaR=52.9,MR=895.7,M1=139.57,M2=493.677;
	double dEdq,E,E1,E2,Gamma,qR,dGammadq,tandelta,dtandeltadq,denom;
	//phaseshifts from P. Estabrooks etal, NPB133, p. 490 (1978).
	
	if(ell==0){
		if(twoI==1){
			a=0.00239;
			reff=-0.00176;
		}
		else{
			a=-0.001;
			reff=-0.00176;  // yes, they are the same values of reff - this is not a typo
		}
		
		denom=((1.0/a)+0.5*reff*q*q);
		tandelta=q/denom;
		dtandeltadq=(tandelta/q)-tandelta*reff*q/denom;
		*delta=atan(tandelta);
		*ddeltadq=dtandeltadq/(1.0+tandelta*tandelta);
	}
	else if(ell==1){
		if(twoI==1){
			E1=sqrt(M1*M1+q*q);
			E2=sqrt(M2*M2+q*q);
			E=E1+E2;
			dEdq=(q/E1)+(q/E2);
			qR=sqrt(Misc::triangle(MR,M1,M2));
			Gamma=GammaR*pow(q/qR,3)*((1.0+pow(qR/R,2))/(1.0+pow(q/R,2)));
			dGammadq=(3.0*Gamma/q)-Gamma*2.0*(q/(R*R))/(1.0+pow(q/R,2));
			tandelta=MR*Gamma/(MR*MR-E*E);
			dtandeltadq=(tandelta*dGammadq/Gamma)+2.0*E*tandelta*dEdq/(MR*MR-E*E);
			*delta=atan(tandelta);
			*ddeltadq=dtandeltadq/(1.0+tandelta*tandelta);
		}
		else{
			*delta=0.0;
			*ddeltadq=0.0;
		}
	}
	else{
		*delta=0.0;
		*ddeltadq=0.0;
	}
	
}

#endif
