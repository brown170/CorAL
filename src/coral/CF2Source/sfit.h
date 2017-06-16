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
#ifndef  __INCLUDE_SFIT_H__
#define  __INCLUDE_SFIT_H__

#include  <cmath>
#include  <cstdio>
#include  <cstring>
#include  <cstdlib>
#include  <string>

#include  "coral.h"

// Li changes
#include  "minimization.h"
using namespace std;

class  CParInfo{
	public :
	bool  fixed;  // 'true' if parameter is not allowed to vary
		// searches are confined to  xmin < x < xmax
	double  xmin,xmax,error,xbar;
	double  bestx,currentx;  // bestx is the value that vave smallest chi^2
	char  *name;
	void  Set(string sname, double  xset, double  error,
		double  xminset, double  xmaxset);
	CParInfo();
	~CParInfo();
};

// Li changes
class  CCF2SFit :  public  CMinimization{
	public :
	/*
	calcflag = 1 if CF is a CCHArray object of specific lx,ly,lz
		and CSourceCalc object makes CCHArray objects
		calcflag = 2 if CF is a C3DArray object 
		and CSourceCalc object makes CCHArray objects
		calcflag =3,4 are used if source functions are calculated through
		intermediate MC lists, but no such implementations yet exist
		calcflag=7, used for comparing 2 CCHArray CF functions
	*/  
		void  SetCalcFlag( int  calcflagset);
	void  SetMCSourceFlag( bool  MCsourceflagset);

	void  SetPar(string parstring, double  value);
	void  SetPar(string parstring, double  value, double  error,double  xmin, double  xmax);
	double GetPar(string parstring);
	void  AddPar(string parstring, double  value, double  error,double  xmin, double  xmax);
	void  PrintPars();
	void  PrintErrorMatrix();
	void  PrintStepMatrix();
	void  FixPar(string parname);
	void  FreePar(string parname);
	void  UseBestPars();
	void  SetL( int  lxset, int  lyset, int  lzset);

// Li changes
	void  ConjugateGradient( int  maxcalls);
	double  fn( double  * x);
	bool  dfn( double  * x);

	void  Metropolis( int  maxcalls);
	void  SteepestDescent( int  maxtries);
	void  Newton( int  maxtries);
	void  UpdateStepMatrix();
	void  InitErrorMatrix();

	// This calculates source functions
	CSourceCalc *sourcecalc;

	// ______________ The objects below depend on calc_flag _____________
		// These are used to store information about the source
	CCHArray *sourceCH;
	C3DArray *source3D;
	CMCList *lista;
	CMCList *listb;

	// Wavefunctions or kernels
	CKernel *kernel;
	CKernelWF *kernelwf;
	CWaveFunction *wf;

	// These are arrays for storing correlation functions and errors
	C3DArray *cexp3D;
	C3DArray *cerror3D;
	C3DArray *ctheory3D;
	CCHArray *cexpCH;
	CCHArray *cerrorCH;
	CCHArray *ctheoryCH;
	// _____________________________________________________________________


	void ResetChiSquared();
	CCF2SFit();
	CCF2SFit(CCHArray *sourceCHset,C3DArray *source3Dset,
		CMCList *listaset,CMCList *listbset,
		CKernel *kernelset,CKernelWF *kernelwfset,
		CWaveFunction *wfset,
		C3DArray *cexp3Dset,C3DArray *cerror3Dset,
		C3DArray *ctheory3Dset,CCHArray *cexpCHset,
		CCHArray *cerrorCHset,CCHArray *ctheoryCHset);
	~CCF2SFit();

	protected :
	int  calcflag;
	static   bool  MCsourceflag;  // If the source has a Monte Carlo nature, 
		// i.e., it fluctuates for a given parameter set, set this to true

	CParInfo **par;
	static   int  nmaxpars;
	int  nfreepars,npars;
	double  **ErrorMatrix;
	double  currentchisquared,bestchisquared;

	int  lx,ly,lz;


	CRandom *randy;
	void  SwitchPars( int  ipara, int  iparb);
	void  SwitchValues( double  *a, double  *b);

	int  ncalls;
	double  **StepMatrix;
	void  Init();
	double  GetChiSquared( double  *x);
	void  CalcErrorMatrixFromCurvature( double  **C);
};

class  CCF2SFit_Blast :  public  CCF2SFit{
	public :
	CCF2SFit_Blast(CSourceCalc *scset,C3DArray *cexpset,
		C3DArray *cerrorset,C3DArray *ctheory3Dset,
		CCHArray *ctheoryset,CCHArray *sourceset,
		CKernel *kernelset);
};

class  CCF2SFit_GX1D :  public  CCF2SFit{
	public :
	CCF2SFit_GX1D(CSourceCalc *scset,CCHArray *cexpset,
		CCHArray *cerrorset,CCHArray *ctheoryset,
		CCHArray *sourceset,CKernel *kernelset);
};

class  CCF2SFit_3DGaussian :  public  CCF2SFit{
	public :
	CCF2SFit_3DGaussian(CSourceCalc *scset,C3DArray *cexpset,
		C3DArray *cerrorset,C3DArray *ctheory3Dset,
		CCHArray *ctheoryset,CCHArray *sourceset,
		CKernel *kernelset);
};

#endif
