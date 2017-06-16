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
#ifndef __INCLUDE_MISC_H
#define __INCLUDE_MISC_H

#include <complex>
#include <cstdlib>
#include <cstring>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_poly.h>
#include "sf.h"
#include "random.h"
#include "constants.h"

using namespace std;

namespace Misc{
	void lorentz(double *u,double *p1,double *p1prime);
// These find tensor/vector in cm frame (same as above)
	void BoostToCM(double *u,double **Pi,double **PiTilde);
	void BoostToCM(double *u,double *p,double *ptilde);

// Given tensor/vector in cm frame, this finds value in lab frame
	void Boost(double *u,double **PiTilde,double **Pi);
	void Boost(double *u,double *ptilde,double *p);


  double cgc(double j1,double m1,double j2,double m2,double j,double m);
  double cgc_edmonds(double j1,double m1,double j2,double m2,double j,double m);
  bool comparestrings(char *s1,char *s2);
  double triangle(double m0,double m1,double m2);
  void outsidelong(double *pa,double *pb,
		   double &qinv,double &qout,double &qside,double &qlong);
  double GetQinv(double *pa,double *pb);
  double GetRapidity(double *pa);
  double GetDely(double *pa,double *pb);
  
  complex<double> cexp(complex<double> z);
  complex<double> ceiphi(double phi);
  complex<double> cpow(complex<double> z,complex<double> a);

  int iround(double x);

  int cgc_delta (int x, int y);
  double cgc_factorial (double n);
  double cgc_fractorial (double n,double m);
  double oldcgc(double j1,double m1,double j2,double m2,double j,double m);

	void Cubic(double a0,double a1,double a2,double a3,
		complex<double>& z1,complex<double>& z2,complex<double>& z3);
  
// Quartic routines provided by Andrew Steiner
	double signswitch(double a, double b);
	void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> &z1, complex<double> &z2,complex<double> &z3,complex<double> &z4);
	void Quartic(double a0,double a1,double a2,double a3,double a4,complex<double> *z);
	// this corresponds to z0=x[0], z1=x[1]+ix[2], z2=x[1]-ix[2]
	void CubicResolvant(double r,double s,double t,double x[],double &d);
	
	int CubicReal(double a0,double a1,double a2,double a3,double *x);
	void CubicComplex(double a0,double a1,double a2,double a3,complex<double> &z1,complex<double> &z2,complex<double> &z3);
  void Pause();
	void Pause(int seconds);
	
	
	
  double CalcDelta_FromSqWells(int ell,double mu,int nwells,double q,double *V0,double *r);
};

#endif
