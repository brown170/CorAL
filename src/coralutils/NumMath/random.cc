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
#ifndef __INCLUDE_CRANDOM_CC
#define __INCLUDE_CRANDOM_CC

#include "random.h"

using namespace std;

CRandom::CRandom(int seed){
  //Choose the algorithm, see gsldoc.ps (page 178)
  //randy=gsl_rng_alloc(gsl_rng_taus); // crappy but fast
  //randy=gsl_rng_alloc(gsl_rng_ranlxd2); // for anal-retentive types
  //randy=gsl_rng_alloc(gsl_rng_ran2); // Not really double precision numbers
  //randy=gsl_rng_alloc(gsl_rng_knuthran2); // sort of slow
  randy=gsl_rng_alloc(gsl_rng_ranlxd1); // Just right

  gsl_rng_set(randy,seed);
}

void CRandom::reset(int seed){
  gsl_rng_set(randy,seed);
}

double CRandom::ran(void){
  return gsl_rng_uniform(randy);
}

long unsigned int CRandom::iran(unsigned long int imax){
  return gsl_rng_uniform_int(randy,imax);
}

double CRandom::gauss(void){
  return gsl_ran_ugaussian(randy);
}

double CRandom::ran_exp(void){
  return -log( ran() );
}

void CRandom::gauss2(double *randy1,double *randy2){
  double x,y,r2,r,c,s;
TRY_AGAIN:
  x=1.0-2.0*gsl_rng_uniform(randy);
  y=1.0-2.0*gsl_rng_uniform(randy);
  r2=x*x+y*y;
  if(r2>1.0) goto TRY_AGAIN;
  r=sqrt(r2);
  c=x/r;
  s=y/r;
  *randy1=c*sqrt(-2.0*log(r2));
  *randy2=(s/c)**randy1;
}

void CRandom::generate_boltzmann(double mass,double T,double *p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3,a,b,c;
  double pmag,ctheta,stheta,phi,pgauss;
  const double NONRELCUTOFF=0.08; // for T/m < cutoff, use non. rel. method
  if(T/mass> NONRELCUTOFF){
  GB_TRYAGAIN:
    r1=ran();
    r2=ran();
    r3=ran();
		a=-log(r1); b=-log(r2); c=-log(r3);
		pmag=T*(a+b+c);
    p[0]=sqrt(pmag*pmag+mass*mass);
    if(ran()>exp((pmag-p[0])/T)) goto GB_TRYAGAIN;
    ctheta=(a-b)/(a+b);
    stheta=sqrt(1.0-ctheta*ctheta);
    phi=T*T*pow(a+b,2)/(pmag*pmag);
    phi=2.0*PI*phi;
    p[3]=pmag*ctheta;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
  }
  else{
    pgauss=sqrt(mass*T);
    gauss2(&p[1],&p[2]);
    p[1]*=pgauss;
    p[2]*=pgauss;
    p[3]=pgauss*gauss();
    p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
  }
}

/*
void CRandom::generate_boltzmann(double mass,double T,double *p){
  const double PI=4.0*atan(1.0);
  double r1,r2,r3;
  double pmag,ctheta,stheta,phi,pgauss;
  const double NONRELCUTOFF=0.08; // for T/m < cutoff, use non. rel. method
  if(T/mass> NONRELCUTOFF){
  GB_TRYAGAIN:
    r1=ran();
    r2=ran();
    r3=ran();
    pmag=-T*log(r1*r2*r3);
    p[0]=sqrt(pmag*pmag+mass*mass);
    if(ran()>exp((pmag-p[0])/T)) goto GB_TRYAGAIN;
    ctheta=log(r1/r2)/log(r1*r2);
    stheta=sqrt(1.0-ctheta*ctheta);
    phi=T*T*pow(log(r1*r2),2)/(pmag*pmag);
    if(phi>1.0){
      printf("phi=%g, out of range\n",phi);
      exit(1);
    }
    phi=2.0*PI*phi;
    p[3]=pmag*ctheta;
    p[1]=pmag*stheta*cos(phi);
    p[2]=pmag*stheta*sin(phi);
  }
  else{
    pgauss=sqrt(mass*T);
    gauss2(&p[1],&p[2]);
    p[1]*=pgauss;
    p[2]*=pgauss;
    p[3]=pgauss*gauss();
    p[0]=sqrt(mass*mass+p[1]*p[1]+p[2]*p[2]+p[3]*p[3]);
  }
}
*/

#endif
