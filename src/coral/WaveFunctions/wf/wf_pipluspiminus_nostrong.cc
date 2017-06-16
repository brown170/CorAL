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
#ifndef __CWAVEFUNCTION_WF_PIPLUSPIMINUS_NOSTRONG_CC__
#define __CWAVEFUNCTION_WF_PIPLUSPIMINUS_NOSTRONG_CC__
#include "wavefunction.h"

CWaveFunction_pipluspiminus_nostrong::CWaveFunction_pipluspiminus_nostrong(string  parsfilename){
	int iq,ichannel,*I;
	double q;

	ParsInit(parsfilename);

	m1=MPI;
	m2=MPI;
	IDENTICAL=0;
	q1q2=-1;
	if(COULOMB==false) q1q2=0;
	nchannels=0;

	InitArrays();
	InitWaves();

	printf("pipluspiminus_nostrong wf initialized\n");
}

double CWaveFunction_pipluspiminus_nostrong::CalcPsiSquared(int iq,double r,double ctheta){
	double psisquared;
	double delta_s,delta_p,delta_d;
	complex<double> psi0;

	if(iq>=nqmax){
		psisquared=1.0;
	}
	else{
		psi0=planewave[iq]->planewave(r,ctheta);
		psisquared=real(psi0*conj(psi0));
		//psisquared*=RelativisticCorrection(r,iq);
	}
	return psisquared;
}

#endif
