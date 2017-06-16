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
#ifndef SOUTWOGAUSSIAN_H
#define SOUTWOGAUSSIAN_H

#include "soubase.h"
#include "objects1d.h"

#ifndef SWIG
enum TwoGaussianParameter {Fraction1,Sigma1,Sigma2};
#endif

class C2GaussianSource: public CSourceFtnBase, public CObject1d
{

 public:
  C2GaussianSource();
  ~C2GaussianSource();

  double getValue(double r) const;
  double getError(double r) const{return 0.;}
  
  void SetParameter(TwoGaussianParameter ParameterName,double value);

  // read/write to parameter map
  virtual bool Read(const parameterMap& s);
  virtual bool Write(parameterMap& s);

 private:
  double mFraction1;
  double mSigma1;
  double mSigma2;
  double Function(double r);
};

#endif
