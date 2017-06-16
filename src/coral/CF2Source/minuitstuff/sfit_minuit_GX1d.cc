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
#ifndef __INCLUDE_SFIT_MINUIT_GX1D_CC__
#define __INCLUDE_SFIT_MINUIT_GX1D_CC__

CCF2S_Minuit_GX1D::CCF2S_Minuit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
				 CCHArray *cerrorset,CCHArray *ctheoryset,
				 CCHArray *sourceset,CKernel *kernelset){
  ndim=1;

  // npars is different for different subclasses
  npars=5;
  //
  sourcecalc=scset;
  cexp=cexpset;
  cerror=cerrorset;
  ctheory=ctheoryset;
  source=sourceset;
  kernel=kernelset;
  if(pars!=NULL) delete [] pars;
  pars=new CMNPars[npars];
  if(xval!=NULL) delete [] xval;
  xval=new double[npars];


  // initialization of pars is also unique to given subclass
  pars[0].Set("lambda",parameter::getD(sourcecalc->spars,"lambdaG",0.6),1.0,0,0.0);
  pars[1].Set("Xfrac",parameter::getD(sourcecalc->spars,"Xfrac",0.5),1.0,0,0);
  pars[2].Set("R",parameter::getD(sourcecalc->spars,"R",5),1.0,0.0,20.0);
  pars[3].Set("X",parameter::getD(sourcecalc->spars,"X",10),1.0,0.0,20.0);
  pars[4].Set("a",parameter::getD(sourcecalc->spars,"a",5),1.0,0,0);

  InitMinuit();
  FixPar(4);
}
#endif
