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
#include <cmath>
#include "constants.h"
#include "message.h"
#include "sou1d_gauss.h"

using namespace std;

//source base class

/// read from parameter map
bool CGaussianSource::Read(const parameterMap& s){
    mSigma = parameter::getD(s,"width",0.);
    mAmplitude = parameter::getD(s,"height",0.);
    mOffset = parameter::getD(s,"offset",0.);
    return CSourceFtnBase::Read(s)&&CObject1d::Read(s);
}

/// write to CCommandOptions map
bool CGaussianSource::Write(parameterMap& s){
    CSourceFtnBase::Write(s)&&CObject1d::Write(s);
    parameter::set(s,"width",mSigma);
    parameter::set(s,"height",mAmplitude);
    parameter::set(s,"offset",mOffset);
    return true;
}


void CGaussianSource::SetParameter(GaussianParameter ParameterName,double value)
{
  switch(ParameterName)
    {
    case Amplitude:
      mAmplitude=value;
      break;
    case Sigma:
      mSigma=value;
      break;
    case Offset:
      mOffset=value;
      break;
    default:
      MESSAGE<<"How the #$^$#^ did you get to here?"<<ENDM_FATAL;
    }
}


double CGaussianSource::Function(double r)
{
  double tmp;
  tmp = (r-mOffset)/mSigma;
  //One must make sure that the function is a *two* particle source function.  That is why this is not exactly the standard gaussian.
  //  return mAmplitude/pow(mSigma*mSigma*2*PI,3.0/2.0)*exp(-tmp*tmp/2);
  return mAmplitude/pow(mSigma*mSigma*4.0*PI,3.0/2.0)*exp(-tmp*tmp/4.0);
}
