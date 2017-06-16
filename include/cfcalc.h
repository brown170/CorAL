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
#ifndef __INCLUDE_CFCALC_H__
#define __INCLUDE_CFCALC_H__
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

#include "arrays.h"

using namespace std;

namespace CFCalc{
  double GetChiSquared(C3DArray *CFexp,C3DArray *Error,C3DArray *CFtheory);
  double GetChiSquared(int lx,int ly,int lz,CCHArray *CFexp,CCHArray *Error,CCHArray *CFtheory);
	double GetChiSquared(CCHArray *CFexp,CCHArray *CFtheory,CCHArray *error);
};

#endif
