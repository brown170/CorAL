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
#include "integratevec.h"
#include "message.h"
#include "cfortran.h" 
#include "utils.h"
#include <cmath>
#include <iostream>
using namespace std;

#define nint iround<double>

//----------------------------------------------
// Fortran wrapper stuff
//----------------------------------------------
PROTOCCALLSFSUB17(DCUHRE,dcuhre,INT,INT,DOUBLEV,DOUBLEV,INT,INT,ROUTINE,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,DOUBLEV,PINT,PINT,DOUBLEV)
#define DCUHRE(NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,WORK) CCALLSFSUB17(DCUHRE,dcuhre,INT,INT,DOUBLEV,DOUBLEV,INT,INT,ROUTINE,DOUBLE,DOUBLE,INT,INT,INT,DOUBLEV,DOUBLEV,PINT,PINT,DOUBLEV,NDIM,NUMFUN,A,B,MINPTS,MAXPTS,FUNSUB,EPSABS,EPSREL,KEY,NW,RESTAR,RESULT,ABSERR,NEVAL,IFAIL,WORK)

//----------------------------------------------
// Global pointers to classes and member functions  
// of that class, so can integrate using member
// functions, solving the callback issue.  Also the 
// wrapper function that wraps the calls
// to the member function.  We'll pass this to the 
// integrator if needed.
//----------------------------------------------
void* Global_CIntegrateVector_ClassPtr;
void (*Global_CIntegrateVector_FuncPtr)(void*,int*,double*,int*,double*);
void Global_CIntegrateVector_FuncWrapper(int* i, double* x, int* j, double* y){
    return Global_CIntegrateVector_FuncPtr(Global_CIntegrateVector_ClassPtr,i,x,j,y);
}

//----------------------------------------------
// default constructor
//----------------------------------------------
CIntegrateVector::CIntegrateVector(void){
    _ndim=0;
    _numfunc=0;
    _lowerlimits=NULL;
    _upperlimits=NULL;
    _minpts=10;
    _maxpts=100;
    _abserr=1e-14;
    _relerr=1e-6;
    _key=0;
    _restart=0;
    _results=NULL;
    _error=NULL;
    _work=NULL;      
    _numftncalls=127;
}

//----------------------------------------------
// default destructor -- deletes all arrays
//----------------------------------------------
CIntegrateVector::~CIntegrateVector(void){
    DeleteWorkArray();
    DeleteResults();
    DeleteLimits();
}

//----------------------------------------------
// check _ndim
//----------------------------------------------
void CIntegrateVector::CheckNDim(void){
    if ((1>_ndim)||(_ndim>15)) {MESSAGE<<"Bad NDIM"<<ENDM_FATAL; exit(-1);}
}

//----------------------------------------------
// check _key
//----------------------------------------------
void CIntegrateVector::CheckKey(void){
    if ((0>_key)||(_key>4)) {MESSAGE<<"Bad KEY"<<ENDM_FATAL; exit(-1);}
}

//----------------------------------------------
// check _maxpts
//----------------------------------------------
void CIntegrateVector::CheckMaxPts(void){
    if (_maxpts < _minpts) {MESSAGE<<"Bad MAXPTS"<<ENDM_FATAL; exit(-1);}
    if (_maxpts<3*_numftncalls) { 
        _maxpts=3*_numftncalls+1;
        MESSAGE << "MAXPTS to low in CIntegrateVector.  Resetting value to "<<
            _maxpts<<ENDM_WARN;
    }
}

//----------------------------------------------
// Creates 1D arrays for the integration limits
//----------------------------------------------
void CIntegrateVector::NewLimits(void){
    if (_lowerlimits != NULL) DeleteLimits();
    if (_upperlimits != NULL) DeleteLimits();
    if (_ndim<=0) {MESSAGE<<"Bad NDIM"<<ENDM_FATAL; exit(-1);}
    _lowerlimits = new double[_ndim];
    _upperlimits = new double[_ndim];
}
//----------------------------------------------
// Deletes 1D arrays of the integration limits
//----------------------------------------------
void CIntegrateVector::DeleteLimits(void){
    if (_lowerlimits != NULL) delete [] _lowerlimits;
    if (_upperlimits != NULL) delete [] _upperlimits;
    _ndim=0;
    _lowerlimits=NULL;
    _upperlimits=NULL;
}

//----------------------------------------------
// Creates 1D arrays for the results vectors
//----------------------------------------------
void CIntegrateVector::NewResults(void){
    if (_results != NULL) DeleteResults();
    if (_numfunc<=0) {MESSAGE<<"Bad NUMFUNC"<<ENDM_FATAL; exit(-1);}
    _results = new double[_numfunc];
    _error = new double[_numfunc];
}

//----------------------------------------------
// Deletes 1D arrays of the results vectors
//----------------------------------------------
void CIntegrateVector::DeleteResults(void){
    if (_results != NULL) delete [] _results;
    if (_error != NULL) delete [] _error;
    _numfunc=0;
    _results=NULL;
    _error=NULL;
}

//----------------------------------------------
// Creates 1D work array
//----------------------------------------------
void CIntegrateVector::NewWorkArray(void){
    if (_work != NULL) DeleteWorkArray();
    SetNW();
    if (_nw<=0) {MESSAGE<<"Bad NW"<<ENDM_FATAL; exit(-1);}
    _work = new double[_nw];
}

//----------------------------------------------
// Deletes 1D work array
//----------------------------------------------
void CIntegrateVector::DeleteWorkArray(void){
    if (_work == NULL) return;
    delete [] _work;
    _nw=0;
    _work=NULL;
}

//----------------------------------------------
// sets _numftncalls
//----------------------------------------------
void CIntegrateVector::SetNumFtnCalls(void){
    if (_ndim == 2) {_numftncalls=65;}
    else if (_ndim == 3) {_numftncalls=127;}
    else {
        _numftncalls = 1 + 4*2*_ndim + 2*_ndim*(_ndim-1) + 4*_ndim*(_ndim-1) +
            4*_ndim*(_ndim-1)*(_ndim-2)/3 + nint(pow(2.,_ndim));
    }
    if (_key == 4) {
        _numftncalls = 1 + 3*2*_ndim + 2*_ndim*(_ndim-1) 
            + nint(pow(2.,_ndim));
    }
}
//----------------------------------------------
// sets _nw
//----------------------------------------------
void CIntegrateVector::SetNW(void){
    int MDIV=1;
    int MAXSUB=(_maxpts-_numftncalls)/(2*_numftncalls)+1;
    _nw=MAXSUB*(2*_ndim+2*_numfunc+2) + 17*_numfunc*MDIV + 1; 
}

//----------------------------------------------
// the main routine(s)
//----------------------------------------------
void CIntegrateVector::Compute(void (*_func)(int*,double*,int*,double*)){
    SetNumFtnCalls();
    CheckMaxPts();
    NewWorkArray();
    DCUHRE(_ndim,_numfunc,_lowerlimits,_upperlimits,_minpts,_maxpts,_func,
        _abserr,_relerr,_key,_nw,_restart,_results,_error,_neval,_ifail,_work);
}
void CIntegrateVector::Compute(void* classptr, void (*_func)(void*,int*,double*,int*,double*)){
    SetNumFtnCalls();
    CheckMaxPts();
    NewWorkArray();
    Global_CIntegrateVector_ClassPtr = classptr;
    Global_CIntegrateVector_FuncPtr = _func;
    DCUHRE(_ndim,_numfunc,_lowerlimits,_upperlimits,_minpts,_maxpts,
        Global_CIntegrateVector_FuncWrapper,
        _abserr,_relerr,_key,_nw,_restart,_results,_error,_neval,_ifail,_work);
}

