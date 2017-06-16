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
#ifndef __INCLUDE_SFIT_GX1D_CC__
#define __INCLUDE_SFIT_GX1D_CC__
#include "sfit.h"

CCF2SFit_GX1D::CCF2SFit_GX1D(CSourceCalc *scset,
			     CCHArray *cexpset,
			     CCHArray *cerrorset,
			     CCHArray *ctheoryset,
			     CCHArray *sourceset,
			     CKernel *kernelset){
  int i,j;

  //
  sourcecalc=scset;
  cexpCH=cexpset;
  cerrorCH=cerrorset;
  ctheoryCH=ctheoryset;
  sourceCH=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass

  AddPar("lambdaG",parameter::getD(sourcecalc->spars,"lambdaG",0.3),
	     0.02,0.0,1.5);
  AddPar("R",parameter::getD(sourcecalc->spars,"R",5),
	     0.2,1.0,12.0);
  AddPar("lambdaX",parameter::getD(sourcecalc->spars,"lambdaX",0.3),
	     0.02,0.0,1.5);
  AddPar("X",parameter::getD(sourcecalc->spars,"X",10.0),
	     0.4,1.0,25.0);
  AddPar("a",parameter::getD(sourcecalc->spars,"a",5.0),
	     0.2,1.0,20.0);

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

