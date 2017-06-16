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
#ifndef __INTEGRATEVEC_H
#define __INTEGRATEVEC_H

// Needed so that can call member functions
enum CallMethod{a_function,a_class};

//---------------------------------------
// Integrator that integrates a vector
// of functions simultaneously.
// It wraps the fortran routine DCHURE.
//---------------------------------------
class CIntegrateVector{

    public:
        // Constructors and Destructors
        CIntegrateVector(void);
        ~CIntegrateVector(void);

        // Initialization
        void SetNDim(int n){_ndim=n;CheckNDim();NewLimits();}
        void SetNumFunc(int n){_numfunc=n;NewResults();}
        void SetMinPts(int n){_minpts=n;}
        void SetMaxPts(int n){_maxpts=n;}
        void SetKey(int n){_key=n;CheckKey();}
        void SetAbsErr(double n){_abserr=n;}
        void SetRelErr(double n){_relerr=n;}
        void SetLimits(int n, double lolim, double uplim)
            {_lowerlimits[n]=lolim;_upperlimits[n]=uplim;}

        // Member access
        int GetNDim(void){return _ndim;}
        int GetNumFunc(void){return _numfunc;}
        int GetMinPts(void){return _minpts;}
        int GetMaxPts(void){return _maxpts;}
        int GetKey(void){return _key;}
        int GetNW(void){return _nw;}
        int GetRestart(void){return _restart;}
        double GetAbsErr(void){return _abserr;}
        double GetRelErr(void){return _relerr;}
        double GetUpperLimit(int n){return _upperlimits[n];}
        double GetLowerLimit(int n){return _lowerlimits[n];}
        int GetNumEvals(void){return _neval;}
        int GetIFail(void){return _ifail;}
        double GetResults(int n){return _results[n];}
        double GetError(int n){return _error[n];}
        
        // Main Routine
        void Compute(void (*_func)(int*,double*,int*,double*));
        void Compute(void* classptr, void (*_func)(void*, int*,double*,int*,double*));
        
    private:
        // Initialization the user shouldn't do
        void SetNW(void);
        void SetNumFtnCalls(void);
        void SetRestart(int n){_restart=n;}
        
        // checkers for integrator parameters
        void CheckNDim(void);
        void CheckKey(void);
        void CheckMaxPts(void);
        
        // Memory management
        void NewLimits(void);
        void DeleteLimits(void);
        void NewResults(void);
        void DeleteResults(void);
        void NewWorkArray(void);
        void DeleteWorkArray(void);

        // Member data that the user supplies somehow
        int _ndim;
        int _numfunc;
        double *_lowerlimits;
        double *_upperlimits;
        int _minpts;
        int _maxpts;
        double _abserr;
        double _relerr;
        int _key;
        int _nw;
        int _restart;

        // Member data we compute ourselves 
        double *_results;
        double *_error;
        int _neval;
        int _ifail;
        double *_work; 
        int _numftncalls;     
};
#endif
