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
#if !defined CHEBYSHEV_H
#define      CHEBYSHEV_H

#if !defined PI
#define PI 3.141592653589793
#endif

#include <iostream>
#include "tnt_array1d.h" /* TNT vector */
#include "tnt_array2d.h" /* TNT matrix */

using namespace TNT;
using namespace std;

double ChebyshevPolynomial(int n, double x);

class CChebyshevApprox1D{

      /* I/O */
      friend ostream& operator<<(ostream &s, const CChebyshevApprox1D &A);
      friend istream& operator>>(istream &s, CChebyshevApprox1D &A);

   public:

      /* Constructors */
      CChebyshevApprox1D(void) {nc = 0;uplim = 0.0;lolim = 0.0;this->RunLoud();};
      CChebyshevApprox1D(
         double upperlimit, double lowerlimit,  
         int ncoeff, double (*func)(double, void*), void *pars)
         {this->Setup(upperlimit,lowerlimit,ncoeff,func,pars);this->RunLoud();};
      CChebyshevApprox1D(
         double upperlimit, double lowerlimit, 
         int ncoeff, const Array1D<double> &controlpts, 
         const Array1D<double> &ftnvals){this->SetupFit(upperlimit, lowerlimit,
           ncoeff, controlpts, ftnvals);this->RunLoud();};
      /* Destructor */
      ~CChebyshevApprox1D(void) {};
      /* Copiers */
      CChebyshevApprox1D Copy(const CChebyshevApprox1D &A);
      CChebyshevApprox1D operator=(const CChebyshevApprox1D &A)
         {this->Copy(A); return *this;};
      /* Member access functions */
      int Ncoeffs(void) const {return nc;};
      double Upperlimit(void) const {return uplim;};
      double Lowerlimit(void) const {return lolim;};
      double Coeff(int i) const {return c[i];};
      Array1D<double> CoeffVector(void) const {return c;};
      /* Initializers */
      void SetParameters(double upperlimit, double lowerlimit, int ncoeff);
      void Setup(
         double upperlimit, double lowerlimit, 
         int ncoeff, double (*func)(double, void*), void *pars);
      void Setup(
         double upperlimit, double lowerlimit, 
         int ncoeff, Array1D<double> coeffs);
      void SetupFit(
         double upperlimit, double lowerlimit, 
         int ncoeff, const Array1D<double> &controlpts, 
         const Array1D<double> &ftnvals);
      /* Colocation Points used to set the approximation */
      double ColocPoint(int i) const;
      /* Approximate value of the function */
      double Val(double x, int m) const;
      double Val(double x) const {return this->Val(x,this->nc);};
      double operator()(double x, int m) const {return this->Val(x,m);};
      double operator()(double x) const {return this->Val(x,this->nc);};
      /* Derivative and integral of the approximation */
      CChebyshevApprox1D Derivative(void) const;
      CChebyshevApprox1D Integral(void) const;
      /* Control amount of output to user */
      void RunQuiet(void) {quiet=true;}
      void RunLoud(void) {quiet=false;}

   protected:

      int nc;
      double uplim, lolim;
      Array1D<double> c;
      bool quiet;

};


class CChebyshevApprox2D{

      /* I/O */
      friend ostream& operator<<(ostream &s, const CChebyshevApprox2D &A);
      friend istream& operator>>(istream &s, CChebyshevApprox2D &A);

   public:

      /* Constructors */
      CChebyshevApprox2D(void);
      CChebyshevApprox2D(
         double upperlimit_x, double lowerlimit_x,  
         double upperlimit_y, double lowerlimit_y,  
         int ncoeff_x, int ncoeff_y, 
         double (*func)(double, double, void*), void *pars)
         {this->Setup(
            upperlimit_x,lowerlimit_x,
            upperlimit_y,lowerlimit_y,
            ncoeff_x,ncoeff_y,func,pars);
         };
      CChebyshevApprox2D(
         double upperlimit_x, double lowerlimit_x, 
         double upperlimit_y, double lowerlimit_y, 
         int ncoeff_x, int ncoeff_y,
         const Array1D<double> &x_controlpts, 
         const Array1D<double> &y_controlpts, 
         const Array2D<double> &ftnvals)
         {this->SetupFit(
            upperlimit_x,lowerlimit_x,
            upperlimit_y,lowerlimit_y,
            ncoeff_x,ncoeff_y,
            x_controlpts,y_controlpts,ftnvals);
         };
      /* Destructor */
      ~CChebyshevApprox2D(void) {}; 
      /* Copiers */
      CChebyshevApprox2D Copy(const CChebyshevApprox2D &A);
      CChebyshevApprox2D operator=(const CChebyshevApprox2D &A)
         {this->Copy(A); return *this;};
      /* Member access functions */
      int Ncoeffs_x(void) const {return nc_x;};
      double Upperlimit_x(void) const {return uplim_x;};
      double Lowerlimit_x(void) const {return lolim_x;};
      int Ncoeffs_y(void) const {return nc_y;};
      double Upperlimit_y(void) const {return uplim_y;};
      double Lowerlimit_y(void) const {return lolim_y;};
      double Coeff(int i, int j) const {return c[i][j];};
      Array2D<double> CoeffMatrix(void) const {return c;};
      /* Initializers */
      void SetParameters(
         double upperlimit_x, double lowerlimit_x,  
         double upperlimit_y, double lowerlimit_y,  
         int ncoeff_x, int ncoeff_y);
      void Setup(
         double upperlimit_x, double lowerlimit_x,  
         double upperlimit_y, double lowerlimit_y, 
         int ncoeff_x, int ncoeff_y, 
         double (*func)(double, double, void*), void *pars);
      void SetupFit(
         double upperlimit_x, double lowerlimit_x, 
         double upperlimit_y, double lowerlimit_y, 
         int ncoeff_x, int ncoeff_y,
         const Array1D<double> &x_controlpts, 
         const Array1D<double> &y_controlpts, 
         const Array2D<double> &ftnvals);
      /* Colocation Points used to set the approximation */
      double ColocPoint_x(int i) const;
      double ColocPoint_y(int i) const;
      /* Ways to get the approximate value of the function */
      double Val(double x, double y, int m_x, int m_y) const;
      double Val(double x, double y) const 
         {return this->Val(x,y,this->nc_x,this->nc_y);};
      double operator()(double x, double y, int m_x, int m_y) const 
         {return this->Val(x,y,m_x,m_y);};
      double operator()(double x, double y) const 
         {return this->Val(x,y,this->nc_x,this->nc_y);};

   protected:

      int nc_x,nc_y;
      double uplim_x, lolim_x;
      double uplim_y, lolim_y;
      Array2D<double> c;

};

#endif
