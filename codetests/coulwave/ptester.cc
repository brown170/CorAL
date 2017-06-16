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
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <complex>

using namespace std;

#include "wavefunction.h"
#include "gslbess.h"
#include "constants.h"

//
// This code was superceeded with the gslbesstester.cc code since we have chosen to 
// use the GSL Bessel function 
//
// DAB 11/5/2009
//

int main(){
  double x;
  int ell;
  printf("Enter ell and x: ");
  scanf("%d %lf",&ell,&x);
  
  printf("(%g,%g) =? (%g,%g)\n",real(hankel(ell,x)),imag(hankel(ell,x)),
	 real(Bessel::hn(ell,x)),imag(Bessel::hn(ell,x)));
  return 0;
}


