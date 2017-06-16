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
#ifndef __INCLUDE_SFIT_MINUIT_3DGAUSSIAN_CC__
#define __INCLUDE_SFIT_MINUIT_3DGAUSSIAN_CC__

CCF2S_Minuit_3DGaussian::CCF2S_Minuit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
				       C3DArray *cerrorset,C3DArray *ctheory3Dset,
				       CCHArray *ctheoryset,CCHArray *sourceset,
				       CKernel *kernelset){
  ndim=3;

  // npars is different for different subclasses
  npars=5;
  //
  sourcecalc=scset;
  cexp3D=cexpset;
  cerror3D=cerrorset;
  ctheory=ctheoryset;
  ctheory3D=ctheory3Dset;
  source=sourceset;
  kernel=kernelset;
  if(pars!=NULL) delete [] pars;
  pars=new CMNPars[npars];
  if(xval!=NULL) delete [] xval;
  xval=new double[npars];

  // initialization of pars is also unique to given subclass
  pars[0].Set("lambda",parameter::getD(sourcecalc->spars,"lambda",1),1.0,0,0);
  pars[1].Set("Rx",parameter::getD(sourcecalc->spars,"Rx",5),1.0,0.0,20.0);
  pars[2].Set("Ry",parameter::getD(sourcecalc->spars,"Ry",5),1.0,0.0,20.0);
  pars[3].Set("Rz",parameter::getD(sourcecalc->spars,"Rz",5),1.0,0.0,20.0);
  pars[4].Set("Xoff",parameter::getD(sourcecalc->spars,"Xoff",0),1.0,0,0);
  //pars[5].Set("Yoff",0.0,1.0,0,0);
  //pars[6].Set("Zoff",0.0,1.0,0,0);

  InitMinuit();
  //FixPar(6);
  //FixPar(5);
  //FixPar(4);
}
#endif

