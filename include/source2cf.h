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
#ifndef __INCLUDE_S2C_H
#define __INCLUDE_S2C_H
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include "random.h"
#include "parametermap.h"
#include "arrays.h"
#include "wavefunction.h"
#include "kernel.h"
#include "sourcecalc.h"

using namespace std;

namespace S2CF{
  void s2c(C3DArray *s,CWaveFunction *wf,C3DArray *cf);
  void s2c(C3DArray *s,CKernelWF *kernel,C3DArray *cf);
  void s2c(C3DArray *s,CWaveFunction *wf,C3DArray *cf);
  void s2c(CCHArray *s,CKernel *kernel,CCHArray *cf);
  void s2c(int Lx,int Ly,int Lz,CCHArray *s,CKernel *kernel,CCHArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CWaveFunction *wf,C3DArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CKernel *kernel,C3DArray *cf);
  void s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf,int NMC);
  void s2c(CMCList *lista,CMCList *listb,CKernelWF *kernel,C3DArray *cf3d);
	void s2c_gauss(CSourceCalc *sourcecalc,CKernelWF *kernel,C3DArray *cf3d);
	void s2c_bowlersinyukov(CSourceCalc *sourcecalc,CKernel *kernel,C3DArray *cf3d);
	void s2c_bosons(CMCList *list,C3DArray *cf3d);
};

#endif
