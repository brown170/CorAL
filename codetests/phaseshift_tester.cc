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

#include "wavefunction.h"
#include "constants.h"
using namespace std;

int main(){
  int iq,I,ell;
  double q,e,delq=50.0,delta,ddeltadq,olddelta;
  double MPI = PionMass;
  printf("Enter I and ell: ");
  scanf("%d %d",&I,&ell);
  olddelta=0.0;
  for(q=0.5*delq;q<2000;q+=delq){
    e=2*sqrt(q*q+MPI*MPI);
    WaveFunctionRoutines::getphaseshift_pipi(I,ell,q,&delta,&ddeltadq);
    delta*=(180.0/PI);
    ddeltadq*=(180.0/PI);
    printf("%1d %1d %6.1f %6.2f %6.3f %6.3f\n",
	   I,ell,e,delta,ddeltadq,(delta-olddelta)/delq);
    olddelta=delta;
  }

  return 0;
}

