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
#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "tnt_array1d_utils.h"
#include "tnt_array1d.h"
#include "jama_svd.h"

//------------------------------------------
// This file contains all of the linear algebra 
// that TNT and JAMA should have came with
//------------------------------------------

using namespace std;
using namespace TNT;
using namespace JAMA;

#define EPSILON 1.0e-15

/// TNT vector = STL vector
template <class T>
Array1D<T> stl2tntVec(const vector<T>& v){
    int N=v.size();
    Array1D<T> a(N);
    for (int i=0;i<N;i++) a[i]=v[i];
    return a;
}

/// STL vector = TNT vector
template <class T>
vector<T> tnt2stlVec(const Array1D<T>& v){
    int N=v.dim();
    vector<T> a;
    for (int i=0;i<N;i++) a.push_back(v[i]);
    return a;
}

/// STL vector of vectors = TNT matrix
template <class T>
vector< vector<T> > tnt2stlMat(const Array2D<T>& m){
    vector< vector<T> > tmp;
    for (int i=0;i<m.dim1();i++){
        vector<T> row;
        for (int j=0;j<m.dim2();j++) row.push_back(m[i][j]);
        tmp.push_back(row);
    }
    return tmp;
}

/// TNT matrix = STL vector of vectors
template <class T>
Array2D<T> stl2tntMat(const vector< vector<T> >& m){
    int ncol=0;
    for(int i=0;i<static_cast<int>(m.size());++i)
        {if (static_cast<int>(m[i].size())>ncol) ncol=m[i].size();}

    Array2D<T> tmp(m.size(),ncol);
    for (int i=0;i<tmp.dim1();++i){
        for (int j=0;j<tmp.dim2();++j) {
            if (j<static_cast<int>(m[i].size())) {tmp[i][j]=m[i][j];} else {tmp[i][j]=static_cast<T>(0);}
        }
    }
    return tmp;
}

/// old version of TNT vector = STL vector
template <class T>
void copy(Array1D<T>& A, const vector<T>& v){
    int N=v.size();
    Array1D<T> a(N);
    for (int i=0;i<N;i++) a[i]=v[i];
    A=a;
}

/// old STL vector = TNT vector
template <class T>
void copy(vector<T>& A, const Array1D<T>& v){
    int N=v.dim();
    vector<T> a;
    for (int i=0;i<N;i++) a.push_back(v[i]);
    A=a;
}

/// old STL vector of vectors = TNT matrix
template <class T>
void copy(vector< vector<T> >& A, const Array2D<T>& m){
    vector< vector<T> > tmp;
    for (int i=0;i<m.dim1();i++){
        vector<T> row(m.dim2());
        for (int j=0;j<m.dim2();j++) row[j]=m[i][j];
        tmp.push_back(row);
    }
    A=tmp;
}

/// TNT matrix = STL vector of vectors 
template <class T>
void copy(Array2D<T>& A, const vector< vector<T> >& m){
    int ncol=0;
    for(int i=0;i<static_cast<int>(m.size());i++){if (static_cast<int>(m[i].size())>ncol) ncol=m[i].size();}

    Array2D<T> tmp(m.size(),ncol);
    for (int i=0;i<tmp.dim1();i++){
        for (int j=0;j<tmp.dim2();j++) {
            if (j<ncol) {tmp[i][j]=m[i][j];} else {tmp[i][j]=static_cast<T>(0);}
        }
    }
    A=tmp;
}


/// Matrix transposed
template <class T>
Array2D<T> transpose(const Array2D<T>& A){
	Array2D<T> B(A.dim2(),A.dim1());
	for (int i=0;i<A.dim1();i++){
		for (int j=0;j<A.dim2();j++){
			B[j][i]=A[i][j];
		}
	}
	return B;
}

/// Matrix * Vector
template <class T>
Array1D<T> operator*(const Array2D<T>& A, const Array1D<T>& v){
	Array1D<T> x(A.dim1());
	for (int i=0;i<A.dim1();i++){
		x[i]=0.0;
		for (int j=0;j<A.dim2();j++){
			x[i]+=A[i][j]*v[j];
		}
	}
	return x;
}

/// Vector * Vector inner product
template <class T>
T inner_prod(const Array1D<T>& A, const Array1D<T>& B){
    T ans=0;
	for (int i=0;i<A.dim1();i++){
        ans=ans+A[i]*B[i];
	}
	return ans;
}

/// T * Vector scalar product
template <class T>
Array1D<T> operator*(const T& x, const Array1D<T>& A){
	Array1D<T> B(A.dim1());
	for (int i=0;i<A.dim1();i++){
        B[i]=x*A[i];
	}
	return B;
}

/// T * Matrix scalar product
template <class T>
Array2D<T> operator*(const T& x, const Array2D<T>& A){
	Array2D<T> B(A.dim1(),A.dim2());
	for (int i=0;i<A.dim1();i++){
        for (int j=0;j<A.dim2();j++){
            B[i][j]=x*A[i][j];
        }
	}
	return B;
}

/// Makes a new Matrix with upper block U and lower block L
template <class T>
Array2D<T> uplow_block(Array2D<T> U, Array2D<T> L){
    if (U.dim2() != L.dim2()) {
        cerr<<"In uplow_block: upper and lower matrices have different dim2!\n";
        exit(-1);
    }
    Array2D<T> M(U.dim1()+L.dim1(),U.dim2());
    for (int j=0;j<U.dim2();j++){
        int i=0;
        for (int k=0;k<U.dim1();k++){
            M[i][j]=U[k][j];i++;
        }
        for (int k=0;k<L.dim1();k++){
            M[i][j]=L[k][j];i++;
        }
    }
    return M;
}

/// Makes a new Matrix with left block L and right block R
template <class T>
Array2D<T> leftright_block(Array2D<T> L, Array2D<T> R){
    if (L.dim1() != R.dim1()) {
        cerr<<"In leftright_block: left and right matrices have different dim1!\n";
        exit(-1);
    }
    Array2D<T> M(L.dim1(),L.dim2()+R.dim2());
    for (int j=0;j<L.dim1();j++){
        int i=0;
        for (int k=0;k<L.dim2();k++){
            M[j][i]=L[j][k];i++;
        }
        for (int k=0;k<R.dim2();k++){
            M[j][i]=R[j][k];i++;
        }
    }
    return M;
}

/// Makes a new Vector with upper block U and lower block L
template <class T>
Array1D<T> uplow_block(Array1D<T> U, Array1D<T> L){
    Array1D<T> V(U.dim1()+L.dim1());
    int i=0;
    for (int k=0;k<U.dim1();k++){
        V[i]=U[k];i++;
    }
    for (int k=0;k<L.dim1();k++){
        V[i]=L[k];i++;
    }
    return V;
}

/// Matrix inversion using the SVD
template <class T>
Array2D<T> svdinvert(const Array2D<T>& A){
	Array2D<T> v,u;
	Array1D<T> w;

	// get svd decomposition of A
	SVD<T> decomp(A);
	decomp.getU(u);
	decomp.getV(v);
	decomp.getSingularValues(w);


    int r=0;
	// prune off singular values that are too small, if needed
	T cutoff_condition=1.0/(EPSILON*A.dim1());
	if (decomp.cond()>cutoff_condition) {
	    cout << "svdinvert: Condition number = "<<decomp.cond();
		cout << ",  failed condition number test, pruning singular values\n";
        cout << "svdinvert: Singular values are:  ";
        for (int i=0;i<w.dim();i++) {cout << w[i] << "  ";}cout <<endl;
		T minw=w[0]/cutoff_condition; //singular values ordered in jama_svd
		for (int i=w.dim()-1;i >=0; i--){
			if (abs(w[i]) < minw) {w[i]=0.0;r++;}
		}
	}
    if (r!=0) cout << "svdinvert: Effective rank = "<<w.dim()-r<<endl;

	// get a square matrix who diagonal values are 1/singular values
	Array2D<T> winv(w.dim(),w.dim());
	winv=0.0;
	for (int i=0;i<w.dim();i++) {if (w[i]!=0.0) {winv[i][i]=1.0/w[i];}}

	// the actual inversion
	return matmult(v,matmult(winv,transpose(u)));
}

//-------------------------------------
// Functions to simplify computation of 
// covariance matrices and error vectors
//-------------------------------------
template <class T>
Array1D<T> CovmtxToErrors(const Array2D<T>& M, const Array2D<T>& cov){
    if ((M.dim2() != cov.dim1())||(M.dim2() != cov.dim2())){
        cerr << "Incompatible dimensions in CovmtxToErrors!"<<endl;
        exit(-1);
    }
        
    Array1D<T> ans(M.dim1(),0.0);
        
    for (int i=0;i<M.dim1();i++){
        ans[i]=0.0;
        for (int j=0;j<M.dim2();j++){
            for (int k=0;k<M.dim2();k++){
                ans[i]+=M[i][j]*M[i][k]*cov[j][k];
            }
        }
        ans[i]=sqrt(ans[i]);
    }
    return ans;
}

template <class T>
Array2D<T> CovmtxToCovMtx(const Array2D<T>& M, const Array2D<T>& cov){
    if ((M.dim2() != cov.dim1())||(M.dim2() != cov.dim2())){
        cerr << "Incompatible dimensions in CovmtxToCovMtx!"<<endl;
        exit(-1);
    }
        
    Array2D<T> ans(M.dim1(),M.dim1(),0.0);
        
    for (int i=0;i<M.dim1();i++){
        for (int j=0;j<M.dim2();j++){
            ans[i][j]=0.0;
            for (int k=0;k<M.dim2();k++){
                for (int l=0;k<M.dim2();k++){
                    ans[i][j]+=M[i][k]*M[j][l]*cov[k][l];
                }
            }
        }
    }
    return ans;
}

template <class T>
Array2D<T> ErrorsToCovmtx(const Array2D<T>& M, const Array1D<T>& err){
    if (M.dim2() != err.dim()){
        cerr << "Incompatible dimensions in ErrorsToCovmtx!"<<endl;
        exit(-1);
    }
    Array2D<T> ans(M.dim1(),M.dim1(),0.0);

    for (int i=0;i<M.dim1();i++){
        for (int j=0;j<M.dim1();j++){
            ans[i][j]=0.0;
            for (int k=0;k<err.dim();k++){
                ans[i][j]+=M[i][k]*M[j][k]*err[k]*err[k];
            }
        }
    }
    return ans;
}

template <class T>
Array2D<T> makeDiag(const Array1D<T>& d){
    Array2D<T> ans(d.dim(),d.dim(),0.);
    for (int i=0;i<d.dim();i++){
        ans[i][i]=d[i];
    }
    return ans;
}

template <class T>
Array2D<T> makeDiagSquared(const Array1D<T>& d){
    Array2D<T> ans(d.dim(),d.dim(),0.);
    for (int i=0;i<d.dim();i++){
        ans[i][i]=d[i]*d[i];
    }
    return ans;
}

#endif
