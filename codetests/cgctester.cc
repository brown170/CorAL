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
#include <complex>

#include "misc.h"
#include "random.h"

using namespace std;
//using namespace Misc;

int main(){
  CRandom randy(12345);
  double J,j1,j2,m1,m2,M;
  int i,n;
  double ans,oldans;
  
  printf("Enter NMC : ");
  scanf("%d",&n);
  //n=50;
  for(i=0;i<n;i++){
    j1=randy.iran(31);
    j2=randy.iran(31);
    j1*=0.5;    j2*=0.5;
    m1=j1-randy.iran(int(2*j1+1));
    m2=j2-randy.iran(int(2*j2+1));
    M=m1+m2;
    J=fabs(j1-j2);
    if(j1>j2){
      J+=randy.iran(int(j2+1));
    }
    else  J+=randy.iran(int(j1+1));
    if(J>(j1+j2+0.00001) || J<fabs(j1-j2)-0.00000001){
      printf("OUCH, J=%g, j1=%g, j2=%g\n",J,j1,j2);
      exit(1);
    }

    ans=Misc::cgc(j1,m1,j2,m2,J,M);
    oldans=0.0; //oldcgc(j1,m1,j2,m2,J,M);
    if(fabs(ans-oldans)>1.0E-5){
      printf("J=%g, j1=%g, m1=%g, j2=%g, m2=%g,         %g =? %g\n",
	     J,j1,m1,j2,m2,ans,oldans);
      exit(1);
    }

    //printf("J=%g, j1=%g, m1=%g, j2=%g, m2=%g,         %g =? %g\n",
    //	   J,j1,m1,j2,m2,ans,oldans);

  }
  return 0;
  
}

