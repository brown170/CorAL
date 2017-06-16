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
/*! \class CRombergIntegrator
\brief A Romberg integrator

This is an implementation of the Romberg integrator from the book numerical
 recipes.  It will integrate functions of type double with a single double
 argument.  The default integration parameters are in the header file, but
they can easily be changed with the alternate constructor.
*/

#include "crombergintegrator.h"
#include <cmath>
#include <iostream>

using namespace std;

/*! The default constructor is just here to make sure that the destructor doesn't destroy memory with a delete of the internal arrays.  This is done by initalizing the arrays to null.
 */
CRombergIntegrator::CRombergIntegrator()
{
  s=0;
  h=0;
}

/*!  This is the constructor to use if you just want to integrate a function, and trust that the default parameters of the Romberg integrator are good enough for you.  In most cases this is a safe bet.  If you plan to torure the integrator, you better learn how to use the constructor with varaible parameters.
 */
CRombergIntegrator::CRombergIntegrator(double(*fun)(double))
{
  mCallMethod=RI_a_function;
  function_fun=fun;
  jmax=JMAX;
  k=K;
  eps=EPS;
  s=new double[jmax];
  h=new double[jmax+1];
}

/*!  This constructor allows you to fully tweek the integrator.  You can only tweek the parameters at the construction of the object.  I didn't feel the need to allow dynamic changing of the parameters.  You can simply make another object with different paramerts.  If this is a big problem for you let me know, and I will make it more fancy.
\param a_eps The fractional accruacy desired.  Be carefull, if you make this too small you will run into machine accuracy.
\param a_jmax The maximum number of divisions the integration will make.  For each increase of 1 in this parameter, the number of division in the function will double.
\param a_k The order of the polynomial that is used in the estimation of the error in the integration.
*/
CRombergIntegrator::CRombergIntegrator(double a_eps, int a_jmax, int a_k,double(*fun)(double))
{
  mCallMethod=RI_a_function;
  function_fun=fun;
  jmax=a_jmax;
  k=a_k;
  eps=a_eps;
  s=new double[jmax];
  h=new double[jmax+1];
}

/*! These next two functions solve the problem when you want to call this integrator from a non static member function.  The problem is that if you send a pointer to the member function that you want to integrate, the compiler will complain since the "this" pointer has to be sent also.  The extra argument allows for this.  You also have to set up a static member in your calling class to pass to this function.  That way this class doesn't have to know anything about the class that is using it.  If you want more information on why this is nessasary, the problem is called "callback" and you can look it up in any good c++ book, or look on the web.
 */
CRombergIntegrator::CRombergIntegrator(void* class_ptr,double(*fun)(void*,double))
{
  mCallMethod=RI_a_class;
  function_class=fun;
  callingClass=class_ptr;
  jmax=JMAX;
  k=K;
  eps=EPS;
  s=new double[jmax];
  h=new double[jmax+1];
}

CRombergIntegrator::CRombergIntegrator(double a_eps, int a_jmax, int a_k,void* class_ptr,double(*fun)(void*,double))
{
  mCallMethod=RI_a_class;
  function_class=fun;
  callingClass=class_ptr;
  jmax=a_jmax;
  k=a_k;
  eps=a_eps;
  s=new double[jmax];
  h=new double[jmax+1];
}

double CRombergIntegrator::function(double x)
{
  switch(mCallMethod)
    {
    case RI_a_function:
      return function_fun(x);
      break;
    case RI_a_class:
      return function_class(callingClass,x);
      break;
    default:
      cerr<<"you should never see this message"<<endl;
      return 0.0;
    }
}

/*! Kill it and all of its memory!!*/
CRombergIntegrator::~CRombergIntegrator()
{
  delete [] s;
  delete [] h;
}

/*! Once you create the object, you have to tell it to compute the integral.  This function will compute the integral and return the result.  You can call this multiple times for a given function without reseting anything.
\param a The start (lower limit) of the integral
\param b The end (upper limit) of the integral
\return The value of the integral over the range of a to b.
*/
double CRombergIntegrator::compute(double a, double b)
{
  double ss,dss;
  int j;
  h[1]=1.0;
  for (j=1;j<=jmax;j++) {
    s[j]=trapzd(a,b,j);
    if (j >= k) {
      polint(&h[j-k],&s[j-k],k,0.0,&ss,&dss);
      if (fabs(dss) <= eps*fabs(ss)) return ss;
    }
    h[j+1]=0.25*h[j];
  }
    cerr<<"Too many steps in routine qromb"<<endl;
    return ss;                    
}

void CRombergIntegrator::polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
  //Given arrays xa[1..n] and ya[1..n], and given a value x, this routine   
  // returns a value y, and an error estimate dy. If P(x) is the polynomial
  // of degree N-1 such that P(xai) =yai ; i=1,...,n, then the returned value
  // y = P(x).                                            
                      
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;
  dif=fabs(x-xa[1]);
  double *c = new double[n+1];
  double *d = new double[n+1];

  for (i=1;i<=n;i++) {  //Here we  nd the index ns of the closest table entry, 
    if ( (dift=fabs(x-xa[i])) < dif)
      {
        ns=i;
        dif=dift;
      }
    c[i]=ya[i];// and initialize the tableau of c's and d's.               
    d[i]=ya[i];
  }
  *y=ya[ns--]; //This is the initial approximation to y.    
  for (m=1;m<n;m++) { //For each column of the tableau,  
    for (i=1;i<=n-m;i++) { //we loop over the current c's and d's and update th\em.   
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) cerr<<"Error in routine polint"<<endl;
      // This error can occur only if two input xa's are (to within roundo ) id\entical.              
      den=w/den;
      d[i]=hp*den; //Here the c's and d's are updated.   
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
  delete [] c;
  delete [] d;
}


double CRombergIntegrator::trapzd(double a, double b, int n)
{
  // This routine computes the nth stage of re nement of an extended trapezoida\l rule. func is input as a pointer to the function to be integrated between lim\its a and b, also input. When called with n=1, the routine returns the crudest \estimate of R b a f(x)dx. Subsequent calls with n=2,3,... (in that sequential o\rder) will improve the accuracy by adding 2 n-2 additional interior points.     
  double x,tnm,sum,del;
  static double s;
  int it,j;
  if (n == 1) {
    return (s=0.5*(b-a)*(function(a)+function(b)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm; // This is the spacing of the points to be added.
    x=a+0.5*del;
    for (sum=0.0,j=1;j<=it;j++,x+=del) sum += function(x);
    s=0.5*(s+(b-a)*sum/tnm); //This replaces s by its re ned value.    
    return s;
  }
}
