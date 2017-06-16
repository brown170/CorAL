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
#include "soubase.h"
#include "sou1d_2gauss.h"
#include <iostream>
#include <cmath>
#include "constants.h"
#include "message.h"

using namespace std;

//source base class

/// read from parameter map
bool C2GaussianSource::Read(const parameterMap& s){
    mSigma1 = parameter::getD(s,"width1",0);
    mSigma2 = parameter::getD(s,"width2",0);
    mFraction1 = parameter::getD(s,"fraction1",0);
    return CSourceFtnBase::Read(s)&&CObject1d::Read(s);
}

/// write to parameter map
bool C2GaussianSource::Write(parameterMap& s){
  CSourceFtnBase::Write(s)&&CObject1d::Write(s);
  parameter::set(s,"width1",mSigma1);
  parameter::set(s,"width2",mSigma2);
  parameter::set(s,"fraction1",mFraction1);
  return true;
}

C2GaussianSource::C2GaussianSource(){
  mSigma1=0;
  mFraction1=0;
  mSigma2=0;
}

C2GaussianSource::~C2GaussianSource()
{
}

double C2GaussianSource::getValue(double r) const
{
  return const_cast<C2GaussianSource*>(this)->Function(r);
}

void C2GaussianSource::SetParameter(TwoGaussianParameter ParameterName,
				    double value){
  switch(ParameterName)
    {
    case Fraction1:
      mFraction1=value;
      break;
    case Sigma1:
      mSigma1=value;
      break;
    case Sigma2:
      mSigma2=value;
      break;
      break;
    default:
      MESSAGE<<"How the #$^$#^ did you get to here?"<<ENDM_FATAL;
    }
}


double C2GaussianSource::Function(double r){
  double tmp,g1,g2,g3;
  //below is the two particle source for a 2 gaussian single particle source
  tmp = r/mSigma1;
  g1= mFraction1*mFraction1/pow(mSigma1*mSigma1*4*PI,3.0/2.0)*exp(-tmp*tmp/4);
  tmp = r/mSigma2;
  g2= (1-mFraction1)*(1-mFraction1)/pow(mSigma2*mSigma2*4*PI,3.0/2.0)*exp(-tmp*tmp/4);
  tmp= r*r/(mSigma1*mSigma1+mSigma2*mSigma2);
  g3= 2*mFraction1*(1-mFraction1)/pow(2*PI*(mSigma1*mSigma1+mSigma2*mSigma2),3.0/2.0)*exp(-tmp/2);
  return g1+g2+g3;
}
