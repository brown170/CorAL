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
#ifndef __INCLUDE_SFIT_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_3DGAUSSIAN_CC__
#include "sfit.h"

class C3DArray;
class CKernel;
class CSourceCalc;

CCF2SFit_3DGaussian::CCF2SFit_3DGaussian(CSourceCalc *scset,
					 C3DArray *cexpset,
					 C3DArray *cerrorset,
					 C3DArray *ctheory3Dset,
					 CCHArray *ctheoryset,
					 CCHArray *sourceset,
					 CKernel *kernelset){
  // npars is different for different subclasses
  //
  sourcecalc=scset;
  cexp3D=cexpset;
  cerror3D=cerrorset;
  ctheoryCH=ctheoryset;
  ctheory3D=ctheory3Dset;
  sourceCH=sourceset;
  kernel=kernelset;
  Init();

  // initialization of pars is also unique to given subclass

  AddPar("lambda",parameter::getD(sourcecalc->spars,"lambda",0.5),
	 0.02,0.0,1.5);
  AddPar("Rx",parameter::getD(sourcecalc->spars,"Rx",5),
	 0.2,1.0,12.0);
  AddPar("Ry",parameter::getD(sourcecalc->spars,"Ry",5),
	 0.2,1.0,12.0);
  AddPar("Rz",parameter::getD(sourcecalc->spars,"Rz",5),
	 0.2,1.0,12.0);
  AddPar("Xoff",parameter::getD(sourcecalc->spars,"Xoff",0.0),
	     0.2,-10.0,10.0);

  InitErrorMatrix();

}
#endif

