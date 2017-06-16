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
#ifndef __GSLMATRIX_H__
#define __GSLMATRIX_H__

#include <cmath>
#include <complex>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_complex.h>

using namespace std;

// REAL MATRICES

class CGSLMatrix_Real{
public:
  int dim;
  // solve y=A*x for x
  void SolveLinearEqs(double *y,double **A,double *x);
  void EigenFind(double **A,double **eigenvec,
		 double *eigenval); // A must be symmetric
  void Invert(double **A,double **Ainv); // A must be symmetric
  CGSLMatrix_Real(int dimset);
  ~CGSLMatrix_Real();
private:
  gsl_eigen_symmv_workspace *w;
  gsl_vector *eval;
  gsl_matrix *evec;
  gsl_vector *g;
  gsl_permutation *p;
  gsl_matrix *m;
  gsl_vector *v;
  double **U;
};

// COMPLEX MATRICES
  
class CGSLMatrix_Complex{
public:
  int dim;
  void SolveLinearEqs(complex<double> *y,complex<double> **A,
		      complex<double> *x);
  void EigenFind(complex<double> **A,complex<double> **eigenvec,
		 double *eigenval); // A must be Hermittian
  void Invert(complex<double> **A,complex<double> **Ainv); // A must  be Herm.
  CGSLMatrix_Complex(int dimset);
  ~CGSLMatrix_Complex();
private:
  gsl_vector *eval;
  gsl_matrix_complex *evec;
  gsl_vector_complex *g;
  gsl_eigen_hermv_workspace *w;
  gsl_permutation *p;
  complex<double> **U;
  gsl_matrix_complex *m;
  gsl_vector_complex *v;
};

#endif
