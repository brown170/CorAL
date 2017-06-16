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
#ifndef SCOTTUTILS_YLM_CC__
#define SCOTTUTILS_YLM_CC__
#include "sf.h"
#include "utils.h"

using namespace std;

double SpherHarmonics::legendre(int ell,double ctheta){
	double answer;
	//printf("check a, ell=%d, ctheta=%g\n",ell,ctheta);
	if(ctheta>1.0) ctheta=1.0;
	if(ctheta<-1.0) ctheta=-1.0;
  return gsl_sf_legendre_Pl(ell,ctheta);
}

complex<double> SpherHarmonics::Ylm(int ell, int m, double theta, double phi){
  double ctheta;
  complex<double> answer;
  complex<double> ci(0.0,1.0);
  ctheta=cos(theta);
  answer=gsl_sf_legendre_sphPlm(ell,abs(m),ctheta)*Misc::ceiphi(static_cast<double>(m)*phi);
  if(m<0) answer *= NEGONE_TO_THE(abs(m)); //pow(-1.0,abs(m));
  return answer;
}

complex<double> SpherHarmonics::Ylm(int ell, int m, double x, double y, double z){
    complex<double> answer;
    double ctheta,phi;
    double r = sqrt(x*x+y*y+z*z);
    if ( r < 1e-10 || fabs(z) < 1e-10 ) ctheta = 0.0;
    else ctheta=z/r;
    phi=atan2(y,x);
    answer=gsl_sf_legendre_sphPlm(ell,abs(m),ctheta)*Misc::ceiphi(static_cast<double>(m)*phi);
    if(m<0) answer *= NEGONE_TO_THE(abs(m));//*pow(-1.0,abs(m));
    return answer;	
}

double SpherHarmonics::ReYlm(int ell, int m, double theta, double phi){
	return real(SpherHarmonics::Ylm(ell,m,theta,phi));
}

double SpherHarmonics::ImYlm(int ell, int m, double theta, double phi){
	return imag(SpherHarmonics::Ylm(ell,m,theta,phi));
}

double SpherHarmonics::ReYlm(int ell, int m, double x,double y,double z){
	return real(SpherHarmonics::Ylm(ell,m,x,y,z));
}

double SpherHarmonics::ImYlm(int ell, int m, double x,double y,double z){
	return imag(SpherHarmonics::Ylm(ell,m,x,y,z));
}

#endif
