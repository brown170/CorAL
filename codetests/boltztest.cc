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
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#include "random.h"

int main(){
  CRandom random(-time(NULL));
  int i,nmc=1000000;
  double mass,pxbar,pybar,pzbar,pbar,ztest,T=100.0;
  double p[4];
  printf("Enter nmc : ");
  scanf("%d",&nmc);
  mass=1.0E-8;
  pxbar=pybar=pzbar=pbar=0.0;
  ztest=0.0;
  for(i=0;i<nmc;i++){
    random.generate_boltzmann(mass,T,p);
    pbar+=sqrt(p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
    pxbar+=p[1];
    pybar+=p[2];
    pzbar+=p[3];
    ztest+=p[1]*p[3];
  }
  pbar=pbar/double(nmc);
  ztest=ztest/double(nmc);
  printf("ztest=%g\n",ztest);
  printf("pbar=%g =? %g\n",pbar,3.0*T);
  pxbar=pxbar/double(nmc);
  pybar=pybar/double(nmc);
  pzbar=pzbar/double(nmc);
  printf("pxbar,pybar,pzbar=(%g,%g,%g) =? 0\n",pxbar,pybar,pzbar);
  return 0;
}
