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
#ifndef __INCLUDE_KERNEL_H
#define __INCLUDE_KERNEL_H

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>
#include <string>
#include "sf.h"
#include "utils.h"
#include "constants.h"
#include "wavefunction.h"
#include "parametermap.h"
using namespace std;

string getKernelFilename( string datadir, int ell, double q );

/**
  *
  *
  */
class CKernel{

 public:
 
  CKernel( string kparsfilename="" );
  virtual ~CKernel();
  virtual double GetValue( int ell, double q, double r ) const;
  double GetValue( int ell, int iq, int ir ) const;
  bool Read( const parameterMap& parameters );
  bool Write( parameterMap& parameters );
  void ReadData( string datadirname );
  void WriteData( string datadirname );
  void Print();
  int GetLMAX();
  double GetDELR();
  double GetDELQ();
  int GetNQMAX();
  int GetNRMAX();
  double GetQ(int iq);
  virtual void Calc( CWaveFunction *wf );
  void Calc_ClassCoul( double ma, double mb, int zazb );
  void Calc_PureHBT();
  bool GetIDENTICAL();
  double GetPsiSquared( int iq, int ir, double ctheta );
  double GetPsiSquared( int iq, double r, double ctheta );
  double GetPsiSquared( double q, double r, double ctheta);

 private:

  bool IDENTICAL;
  int ellmax;
  int nrmax,nqmax;
  double delr,delq;
  double ***kernel;
  double *P;
  void ParsInit( string kparsfilename );
  double CalcPsiSquared( int iq, int ir, double ctheta );
  double CalcPsiSquared( int iq, double r, double ctheta );
  void CalcP( double ctheta );
};

class CKernelExactHBT: public CKernel {
  public:
  
    CKernelExactHBT( string kparsfilename="" ): CKernel( kparsfilename ){}
    double GetValue( int ell, double q, double r ) const 
        {return NEGONE_TO_THE(ell)*Bessel::jn(ell,2.0*q*r/HBARC);}
};


class CKernelWF{

 public:

  double GetPsiSquared(int iq,int ir,int ictheta);
  double GetPsiSquared(int iq,int ir,double ctheta);
  double GetPsiSquared(int iq,double r,double ctheta);
  double GetPsiSquared(double q,double r,double ctheta);
  void Calc(CWaveFunction *wf);
  void ReadData( string datadirname );
  void WriteData( string datadirname );
  double GetDELR();
  double GetDELQ();
  double GetDELCTHETA();
  int GetNQMAX();
  int GetNRMAX();
  int GetNCTHETA();
  bool GetIDENTICAL();
  CKernelWF( string kparsfilename );
  ~CKernelWF();

 private:

  bool IDENTICAL;
  int nctheta,nrmax,nqmax;
  double delr,delq,delctheta;
  double ***wfarray;
  void ParsInit( string kparsfilename );
  
};

#endif
