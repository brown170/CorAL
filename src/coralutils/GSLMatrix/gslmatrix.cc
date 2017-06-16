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
#ifndef __GSL_MATRIX_CC__
#define __GSL_MATRIX_CC__
#include "gslmatrix.h"


CGSLMatrix_Real::CGSLMatrix_Real(int dimset){
	dim=dimset;
	w=NULL;
	eval=NULL;
	evec=NULL;
	g=NULL;
	p=NULL;
	U=NULL;
	m=NULL;
	v=NULL;
}

CGSLMatrix_Real::~CGSLMatrix_Real(){
	if(w!=NULL) gsl_eigen_symmv_free(w);
	if(eval!=NULL) gsl_vector_free(eval);
	if(evec!=NULL) gsl_matrix_free(evec);
	if(g!=NULL) gsl_vector_free(g);
	if(p!=NULL) gsl_permutation_free(p);
	if(m!=NULL) gsl_matrix_free(m);
	if(v!=NULL) gsl_vector_free(v);
	if(U!=NULL){
		for(int i=0;i<dim;i++) delete [] U[i];
		delete [] U;
	}
}

void CGSLMatrix_Real::SolveLinearEqs(double *y,double **A,double *x){
	int i,j,s;
	if(v==NULL) v=gsl_vector_alloc(dim);

	if(m==NULL) m=gsl_matrix_alloc(dim,dim);

	for(i=0;i<dim;i++){
		gsl_vector_set(v,i,y[i]);
		for(j=0;j<dim;j++) gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(g==NULL) g = gsl_vector_alloc (dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);

	gsl_linalg_LU_decomp (m,p,&s);
	gsl_linalg_LU_solve (m, p,v,g);

	for(i=0;i<dim;i++) x[i]=gsl_vector_get(g,i);
}

void CGSLMatrix_Real::EigenFind(double **A,double **eigenvec,
double *eigenval){
	int i,j;
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_symmv_alloc(dim);

	gsl_eigen_symmv(m,eval,evec,w);
	//gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	for(i=0;i<dim;i++){
		eigenval[i]=gsl_vector_get(eval,i);
		for(j=0;j<dim;j++) eigenvec[i][j]=gsl_matrix_get(evec,i,j);
	}

}

void CGSLMatrix_Real::Invert(double **A,double **Ainv){

	int i,j;
	if(m==NULL) m=gsl_matrix_alloc(dim,dim);
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) gsl_matrix_set(m,i,j,A[i][j]);
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_symmv_alloc(dim);

	gsl_eigen_symmv(m,eval,evec,w);
	//gsl_eigen_symmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	if(U==NULL){
		U=new double*[dim];
		for(i=0;i<dim;i++) U[i]=new double[dim];
	}

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) {
			U[i][j]=gsl_matrix_get(evec,j,i);
			Ainv[i][j]=0.0;
		}
	}  
	int k;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++) Ainv[i][j]+=U[k][i]*U[k][j]/gsl_vector_get(eval,k);
		}
	}
}

// HERMITIAN COMPLEX MATRICES


CGSLMatrix_Complex::CGSLMatrix_Complex(int dimset){
	dim=dimset;
	eval=NULL;
	evec=NULL;
	w=NULL;
	g=NULL;
	p=NULL;
	U=NULL;
	m=NULL;
	v=NULL;
}

CGSLMatrix_Complex::~CGSLMatrix_Complex(){
	if(eval!=NULL) gsl_vector_free(eval);
	if(evec!=NULL) gsl_matrix_complex_free(evec);
	if(w!=NULL) gsl_eigen_hermv_free(w);
	if(g!=NULL) gsl_vector_complex_free(g);
	if(p!=NULL) gsl_permutation_free(p);
	if(m!=NULL) gsl_matrix_complex_free(m);
	if(v!=NULL) gsl_vector_complex_free(v);
	if(U!=NULL){
		for(int i=0;i<dim;i++) delete [] U[i];
		delete [] U;
	}
}

void CGSLMatrix_Complex::EigenFind(complex<double> **A,
complex<double> **eigenvec,double *eigenval){
	complex<double> ci(0.0,1.0);
	gsl_complex z;
	//gsl_matrix_complex *m=gsl_matrix_complex_alloc(dim,dim);
	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_complex_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_hermv_alloc(dim);

	gsl_eigen_hermv(m,eval,evec,w);

	//gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	for(i=0;i<dim;i++){
		eigenval[i]=gsl_vector_get(eval,i);
		for(j=0;j<dim;j++){
			z=gsl_matrix_complex_get(evec,i,j);
			eigenvec[i][j]=GSL_REAL(z)+ci*GSL_IMAG(z);
		}
	}

}

void CGSLMatrix_Complex::Invert(complex<double> **A,complex<double> **Ainv){
	complex<double> ci(0.0,1.0);
	gsl_complex z;
	//gsl_matrix_complex *m=gsl_matrix_complex_alloc(dim,dim);
	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	int i,j;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(eval==NULL) eval=gsl_vector_alloc(dim);
	if(evec==NULL) evec=gsl_matrix_complex_alloc(dim,dim);
	if(w==NULL) w=gsl_eigen_hermv_alloc(dim);

	gsl_eigen_hermv(m,eval,evec,w);
	//gsl_eigen_hermv_sort(eval,evec,GSL_EIGEN_SORT_ABS_ASC);

	if(U==NULL){
		U=new complex<double> *[dim];
		for(i=0;i<dim;i++) U[i]=new complex<double>[dim];
	}

	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++) {
			z=gsl_matrix_complex_get(evec,j,i);
			U[i][j]=GSL_REAL(z)+ci*GSL_IMAG(z);
			Ainv[i][j]=0.0;
		}
	}  
	int k;
	for(i=0;i<dim;i++){
		for(j=0;j<dim;j++){
			for(k=0;k<dim;k++)
				Ainv[i][j]+=U[k][i]*conj(U[k][j])/gsl_vector_get(eval,k);
		}
	}

}

void CGSLMatrix_Complex::SolveLinearEqs(complex<double> *y,
complex<double> **A,complex<double> *x){
	complex<double> ci(0.0,1.0);
	int i,j,s;
	gsl_complex z;

	if(m==NULL) m=gsl_matrix_complex_alloc(dim,dim);
	if(v==NULL) v=gsl_vector_complex_alloc(dim);
	for(i=0;i<dim;i++){
		GSL_SET_COMPLEX(&z,real(y[i]),imag(y[i]));
		gsl_vector_complex_set(v,i,z);
		for(j=0;j<dim;j++){
			GSL_SET_COMPLEX(&z,real(A[i][j]),imag(A[i][j]));
			gsl_matrix_complex_set(m,i,j,z);
		}
	}

	if(g==NULL) g = gsl_vector_complex_alloc (dim);
	if(p==NULL) p = gsl_permutation_alloc (dim);


	gsl_linalg_complex_LU_decomp (m, p, &s);
	gsl_linalg_complex_LU_solve (m, p, v, g);

	for(i=0;i<dim;i++){
		z=gsl_vector_complex_get(g,i);
		x[i]=GSL_REAL(z)+ci*GSL_IMAG(z);
	}
}

#endif
