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
#ifndef __INCLUDE_SFIT_BLAST_CC__
#define __INCLUDE_SFIT_BLAST_CC__
#include "sfit.h"

CCF2SFit_Blast::CCF2SFit_Blast(CSourceCalc *scset,C3DArray *cexp3Dset,C3DArray *cerror3Dset, C3DArray *ctheory3Dset,CCHArray *ctheoryset,CCHArray *sourceset,CKernel *kernelset){
	int i,j;
	
	// nfreepars is different for different subclasses
	
	sourcecalc=scset;
	cexp3D=cexp3Dset;
	cerror3D=cerror3Dset;
	ctheoryCH=ctheoryset;
	ctheory3D=ctheory3Dset;
	sourceCH=sourceset;
	kernel=kernelset;
	Init();
	
	// initialization of pars is also unique to given subclass
	
	AddPar("lambda",parameter::getD(sourcecalc->spars,"lambda",1),
		0.01,0.0,1.5);
	AddPar("R",parameter::getD(sourcecalc->spars,"R",13),
		0.05,5.0,20.0);
	AddPar("Tau",parameter::getD(sourcecalc->spars,"Tau",12),
		0.05,5.0,20.0);
	AddPar("DelTau",parameter::getD(sourcecalc->spars,"DelTau",5),
		0.05,0.0,20.0);
	
	
	for(i=0;i<nfreepars;i++){
		for(j=0;j<nfreepars;j++){
			StepMatrix[i][j]=0.0;
			ErrorMatrix[i][j]=0.0;
			if(i==j){
				StepMatrix[i][j]=par[i]->error;
				ErrorMatrix[i][j]=par[i]->error*par[i]->error;
			}
		}
	}
}
#endif

