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
#ifndef __CROMBERG_INTEGRATOR_H
#define __CROMBERG_INTEGRATOR_H

/*!
  \brief A romberg integrator.

  See numerical recipes
*/


//default values -- if you torture the integrator a lot, you may have to vary these with the alternate constructor
#define EPS 1.0e-6  //fractional accuracy desired
#define JMAX 20  // max number of steps
#define K 5  // number of polints used in extrapolation

enum RI_CallMethod{RI_a_function,RI_a_class};

class CRombergIntegrator
{
 public:
  CRombergIntegrator();
  CRombergIntegrator(double(*function)(double));
  CRombergIntegrator(double eps, int jmax, int k, double(*function)(double));
  CRombergIntegrator(void*,double(*function)(void*,double));
  CRombergIntegrator(double eps, int jmax, int k, void*,double(*function)(void*,double));
  ~CRombergIntegrator();

  double compute(double a, double b);

 protected:
  //two methods of getting the function
  double (*function_class)(void*,double);  //pointer for function to be integrated
  void* callingClass;
  double (*function_fun)(double);  //pointer for function to be integrated
  RI_CallMethod mCallMethod;

  double function(double);
  double *s; //array of integral value
  double *h; //array of step sizes 
  int jmax; //max number of steps
  int k; //number of polints used in extrapolation
  double eps;
  void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
  double trapzd(double a, double b, int n);
};

#endif
