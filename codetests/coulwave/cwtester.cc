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
#include "sf.h"
#include "constants.h"

//
// This code needs the function  CWincoming_smallr  which doesn't appear to be 
// declared or implemented anywhere
//
// DAB 11/5/2009
//

void SphericalCW(int ell,double x,double eta,double *FL,double *GL);
void SphericalCWprime(int L,double x,double eta,double *FL,double *GL,
		      double *FLprime,double *GLprime);

int main(){
  double x,eta,abscg,abscg0;
  double FL,GL,rcw,icw;
  int ell;
  printf("Enter ell and x: ");
  scanf("%d %lf",&ell,&x);
  
  for(eta=-3.0;eta<=3.01;eta+=0.2){
    SphericalCW(ell,x,eta,&FL,&GL);
    rcw=real(CWincoming_smallr(ell,x,eta));
    icw=imag(CWincoming_smallr(ell,x,eta));
    printf("eta=%6.3f:  (%g,%g) =? (%g,%g),  %g =? %g, off by %g\n",eta,
	   rcw,icw,FL,GL,sqrt(rcw*rcw+icw*icw),sqrt(FL*FL+GL*GL), 
	   sqrt(rcw*rcw+icw*icw)/sqrt(FL*FL+GL*GL));
  }

  return 0;
}

void SphericalCW(int L,double x,double eta,double *FL,double *GL){
  double expF,expG;
  double *fc,*gc;
  int k=0;
  fc=new double[L+1];
  gc=new double[L+1];

  // This calculates fc and gc arrays for indices L to L+k  
  int dummy=gsl_sf_coulomb_wave_FG_array(0,L,eta,x,fc,gc,&expF,&expG);
  *FL=fc[L]*exp(expF);
  *GL=gc[L]*exp(expG);
  delete [] fc;
  delete [] gc;
}

// This calcules dCW/dx
void SphericalCWprime(int L,double x,double eta,double *FL,double *GL,
		      double *FLprime,double *GLprime){
  double expF,expG;
  double *fc,*gc,*fcp,*gcp;
  int k=0;
  fc=new double[k+1];
  gc=new double[k+1];
  fcp=new double[k+1];
  gcp=new double[k+1];
  // This calculates fc and gc arrays for indices L to L+k  
  int dummy=gsl_sf_coulomb_wave_FGp_array(L,k,eta,x,fc,fcp,gc,gcp,
					 &expF,&expG);
  *FL=fc[0]*exp(expF);
  *GL=gc[0]*exp(expG);
  *FLprime=fcp[0]*exp(expF);
  *GLprime=gcp[0]*exp(expG);
  delete [] fc;
  delete [] gc;
  delete [] fcp;
  delete [] gcp;
}
