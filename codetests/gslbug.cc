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
#include <cstdio>
#include <cmath>
#include <gsl/gsl_sf.h>

int main(){
  int L=2,dummy;
  double x,eta=0.01;
  double expF,expG;
  double dphi0F,dphi1F,dphi2F,dphi0G,dphi1G,dphi2G;
  double F0,F1,F2,G0,G1,G2;
  double *fc,*gc;
  fc=new double[L+1];
  gc=new double[L+1];

  printf("eta=%g\n",eta);
  printf("x ____ F(L=0)  G(L=0) ___  F(L=1) G(L=1) ___ F(L=2) G(L=2) ___\n");
  for(x=0.1;x<7;x+=0.05){
    dummy=gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
    expF=exp(expF);
    expG=exp(expG);
    F0=expF*fc[0];
    F1=expF*fc[1];
    F2=expF*fc[2];
    G0=expG*gc[0];
    G1=expG*gc[1];
    G2=expG*gc[2];
    printf("%4.2f   %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f   %7.3f %7.3f\n",
	            x,  F0,   G0,     dphi0F,dphi0G,F1,G1, dphi1F,dphi1G,F2,G2, dphi2F,dphi2G);
  }
  delete [] fc;
  delete [] gc;

  return 0;
}








