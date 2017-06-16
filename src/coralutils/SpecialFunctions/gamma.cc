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
#include <math.h>
#include "special_func.h"
#include <iostream>
#include <complex>

using namespace std;

/*
   ...gamma(z)...

   Gamma function with complex argument
   translated from FORTRAN by D. Brown 2/11/04
   originally from P. Danielewicz's ucofpak
*/
double_complex Gamma(double_complex z)
{
	const double PI=4.0*atan(1.0);
    int l,k;
    double_complex u,v,h,s,fk,cdgamma;
    static double g[16]={
        41.624436916439068e0,-51.224241022374774e0,
        +11.338755813488977e0,-0.747732687772388e0, 
        +0.008782877493061e0, -0.000001899030264e0,
        +0.000000001946335e0, -0.000000000199345e0,
        +0.000000000008433e0, +0.000000000001486e0, 
        -0.000000000000806e0, +0.000000000000293e0,
        -0.000000000000102e0, +0.000000000000037e0,
        -0.000000000000014e0, +0.000000000000006e0
    };
    double x;

    u=z;
    while(true) {
        x=u.real();
 
        if ( x >= 1.) {
            v=u;
            l=3;
        } else if (x >= 0. ) {
            v=u+1.0;
            l=2;
        } else {
            v=1.0-u;
            l=1;
        }

        h=double_complex(1.0,0.);
        s=g[0];
        for (k = 2; k <= 16; k++) {
        fk=k-2.0;
        h=((v-(fk+1.0))/(v+fk))*h;
        s=s+g[k-1]*h;
        }
        h=v+4.5;
        cdgamma=2.506628274631001*exp((v-0.5)*log(h)-h)*s;

        if (l < 0) {
            return PI/(sin(PI*u)*cdgamma);
        } else if (l == 0) {
            return cdgamma/u;
        } else {
            return cdgamma;
        }


        u=z+1.0;
    }

}

/*
   ...ln(gamma(xx))...

   log of the Gamma function (from Numerical Recipes in C)
*/
double LnGamma(double xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
/*
   ...gamma(xx)...
   
   watch for floating overflows
*/
double Gamma(double xx)
{
  return exp(LnGamma(xx));
}
/*
   ...Factorial(n) = n!...
   
   adapted from Numerical Recipes in C
*/
double Factorial(int n)
{
  static int ntop=4;
  static double answertable[33]={1.0,1.0,2.0,6.0,24.0};
  int i;

  if (n<0) {cerr << "Negative factorial in function Factorial!"<<endl;}
  if (n>32) return Gamma((double)n+1.0);
  while (ntop<n) {
    i=ntop++;
    answertable[ntop]=answertable[i]*ntop;
  }
  return answertable[n];
}
/*
   ... Binomial Coefficient ...
   
   adapted from Numerical Recipes in C
*/
double BinomialCoeff(int n, int k)
{
  return floor(0.5+exp(LnFactorial(n)-LnFactorial(k)-LnFactorial(n-k)));
}
/*
   ... ln(n!)...
   
   adapted from Numerical Recipes in C
*/
double LnFactorial(int n)
{
  static double answertable[101];
 
  if (n < 0) cerr<<"Negative factorial in function LnFactorial";
  if (n <= 1) return 0.0;
  if (n <= 100) return answertable[n] ? answertable[n] : (answertable[n]=LnGamma((double)n+1.0));
  else return LnGamma((double)n+1.0);
}
// Headers for series solution of Incomplete Gamma Functions P and Q
void contfrac_IncGammaQ(double *gammcf, double a, double x, double *gln);
void series_IncGammaP(double *gamser, double a, double x, double *gln);
/*
   ... Incomplete Gamma Function P ...
   
   adapted from Numerical Recipes in C
*/
double IncGammaP(double a, double x)
{
  double gamser,gammcf,gln;
 
  if (x < 0.0 || a <= 0.0) cerr<<"Invalid arguments in function IncGammaP"<<endl;
  if (x < (a+1.0)) {
    series_IncGammaP(&gamser,a,x,&gln);
    return gamser;
  } else {
    contfrac_IncGammaQ(&gammcf,a,x,&gln);
    return 1.0-gammcf;
  }
}
/*
   ... Incomplete Gamma Function Q ...
   
   adapted from Numerical Recipes in C
*/
double IncGammaQ(double a, double x)
{
  double gamser,gammcf,gln;
 
  if (x < 0.0 || a <= 0.0) cerr<<"Invalid arguments in routine IncGammaQ"<<endl;
  if (x < (a+1.0)) {
    series_IncGammaP(&gamser,a,x,&gln);
    return 1.0-gamser;
  } else {
    contfrac_IncGammaQ(&gammcf,a,x,&gln);
    return gammcf;
  }
}
/*
   ... series evaulation of Incomplete Gamma Function P ...
   
   adapted from Numerical Recipes in C
*/
#define INCGAMP_ITMAX 100
#define INCGAMP_EPS 3.0e-7
 
void series_IncGammaP(double *gamser, double a, double x, double *gln)
{
  int n;
  double sum,del,ap;
 
  *gln=LnGamma(a);
  if (x <= 0.0) {
    if (x < 0.0) cerr<<"x less than 0 in routine series_IncGammaP"<<endl;
    *gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=1;n<=INCGAMP_ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*INCGAMP_EPS) {
        *gamser=sum*exp(-x+a*log(x)-(*gln));
        return;
      }
    }
    cerr<<"a too large, INCGAMP_ITMAX too small in routine series_IncGammaP"<<endl;
    return;
  }
}
#undef INCGAMP_ITMAX
#undef INCGAMP_EPS
/*
   ... continued fraction evaulation of Incomplete Gamma Function Q ...
   
   adapted from Numerical Recipes in C
*/
#define INCGAMQ_ITMAX 100
#define INCGAMQ_EPS 3.0e-7
#define INCGAMQ_FPMIN 1.0e-30
 
void contfrac_IncGammaQ(double *gammcf, double a, double x, double *gln)
{
  int i;
  double an,b,c,d,del,h;
 
  *gln=LnGamma(a);
  b=x+1.0-a;
  c=1.0/INCGAMQ_FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=INCGAMQ_ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < INCGAMQ_FPMIN) d=INCGAMQ_FPMIN;
    c=b+an/c;
    if (fabs(c) < INCGAMQ_FPMIN) c=INCGAMQ_FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < INCGAMQ_EPS) break;
  }
  if (i > INCGAMQ_ITMAX) cerr<<"a too large, INCGAMQ_ITMAX too small in contfrac_IncGammaQ"<<endl;
  *gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef INCGAMQ_ITMAX
#undef INCGAMQ_EPS
#undef INCGAMQ_FPMIN
