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
#ifndef __INCLUDE_MINUIT_H__
#define __INCLUDE_MINUIT_H__

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <string>

#include "sfit.h"
#include "minuitfcn.h"
#include "cfortran.h"
#include "cminuit.h"

PROTOCCALLSFSUB3(MNINIT,mninit,PINT,PINT,PINT);
PROTOCCALLSFSUB7(MNPARM,mnparm,INT,STRING,DOUBLE,DOUBLE,PDOUBLE,PDOUBLE,PINT);
PROTOCCALLSFSUB4(MNCOMD,mncomd,ROUTINE,STRING,PINT,ROUTINE);
PROTOCCALLSFSUB8(MNCONT,mncont,ROUTINE,INT,INT,INT,PDOUBLE,PDOUBLE,
		 PINT,ROUTINE);
PROTOCCALLSFSUB6(MNSTAT,mnstat,PDOUBLE,PDOUBLE,PDOUBLE,PINT,PINT,PINT);
PROTOCCALLSFSUB7(MNPOUT,mnpout,INT,PSTRING,PDOUBLE,PDOUBLE,PDOUBLE,
		 PDOUBLE,PINT);
PROTOCCALLSFSUB2(MNFIXP,mnfixp,INT,PINT);
PROTOCCALLSFSUB1(MNFREE,mnfree,INT);

class CMNPars {
public:
  char name[20];
  double value;
  double error;
  double min;
  double max;
  void Set(string parname,double value,double error,
	   double min,double max);
};

class CCF2S_Minuit{
public:
  static int ndim;  // =1 or 3 for 1/3 dimensional CF
  static const int MAXLINE=80;
  static int npars,dummy,Ncalls;
  static double *xval,*grad;
  static CMNPars *pars;
  static CSourceCalc *sourcecalc;
  static CKernel *kernel;


  static CCHArray *ctheory;
  // These are used for 3D Calculations
  static C3DArray *cexp3D;
  static C3DArray *cerror3D;
  static C3DArray *ctheory3D;
  //These are used for 1D Calculations
  static int lx,ly,lz;
  static CCHArray *cexp;
  static CCHArray *cerror;
  static CCHArray *source;

  static void InitMinuit();
  static void CalcChiSquare(int npar,double* grad,double* fcnval,double* xval,
			    int iflag,void* futil);

  static void CalcChiSquare3D(int npar,double* grad,double* fcnval,double* xval,
			    int iflag,void* futil);
  static void CalcChiSquare1D(int npar,double* grad,double* fcnval,double* xval,
			    int iflag,void* futil);
  static void Scan(int ipar,int npts,double start,double end);
  static void Minimize();
  static void Minimize(int maxcalls);
  static void Minimize(int maxcalls,double tolerance);
  static void Migrad();
  static void Migrad(int maxcalls);
  static void Migrad(int maxcalls,double tolerance);
  static void Simplex();
  static void Simplex(int maxcalls);
  static void Simplex(int maxcalls,double tolerance);
  static void Minos();
  static void Minos(int ipar);
  static void StratLevel(int istrategy);
  static void Mnstat();
  static void ViewPars();
  static void SetError(double error);
  static void ErrorMatrix();
  static void Contour(int iparx,int ipary,int npts,
		      double *xcontour,double *ycontour);
  static void SetPar(int ipar,char *name,double value,double error,
		     double min,double max);
  static void FixPar(int ipar);
  static void FreePar(int ipar);
  CCF2S_Minuit();
};


class CCF2S_Minuit_3DGaussian : public CCF2S_Minuit{
public:
  CCF2S_Minuit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
		     C3DArray *cerrorset,C3DArray *ctheory3Dset,
		     CCHArray *ctheoryset,CCHArray *sourceset,
		     CKernel *kernelset);
};

class CCF2S_Minuit_Blast : public CCF2S_Minuit{
public:
  CCF2S_Minuit_Blast(CSourceCalc *scset,C3DArray *cexpset,
		     C3DArray *cerrorset,C3DArray *ctheory3Dset,
		     CCHArray *ctheoryset,CCHArray *sourceset,
		     CKernel *kernelset);
};

class CCF2S_Minuit_GX1D : public CCF2S_Minuit{
 public:
  CCF2S_Minuit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
	       CCHArray *cerrorset,CCHArray *ctheoryset,
	       CCHArray *sourceset,CKernel *kernelset);
};
#endif
