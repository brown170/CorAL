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
#include <iostream>
#include "chebyshev.h"
#include "tnt_array1d.h"  /* TNT vector */
#include "tnt_array2d.h" /* TNT matrix */
#include "jama_svd.h"  /* JAMA SVD */ 
#include "linalg.h" /* extras not in TNT or JAMA */

using namespace TNT;
using namespace JAMA;

//----------------------------------------------------------//
//                                                          //
//  The Chebyshev Polynomial Itself                         //
//                                                          //
//----------------------------------------------------------//
/** 
 *  \fn double ChebyshevPolynomial(int n, double x)
 *  \brief Returns the Chebyshev Polynomial \f$ T_n(x) \f$.
 *  \param n order of the polynomial you want
 *  \param value you want the polynomial evaulated at  
 */
double ChebyshevPolynomial(int n, double x){
   if (n!=0) {
      CChebyshevApprox1D Tk;   
      Array1D<double> coeffs(n+1);
      coeffs=0.0;
      coeffs[n]=1.0;
      Tk.Setup(1.0,-1.0,n+1,coeffs);
      return Tk(x);
   } else {return 1.0;} 
}

//----------------------------------------------------------//
//                                                          //
// Main class definition file for CChebyshevApprox1D.       //
// Documentation is written using Doxygen.                  //
//                                                          //
//----------------------------------------------------------//
/** 
 *  \fn CChebyshevApprox1D::CChebyshevApprox1D(void) 
 *  \brief Plain vanilla constructor for the 1D approximation.  
 */

/**
 *  \fn CChebyshevApprox1D::CChebyshevApprox1D(double upperlimit, double lowerlimit, int ncoeff, double (*func)(double, void*), void *pars)
 *  \brief This routine is both the constructor and setup function bundled 
 *         together.
 *  \param upperlimit upper limit of the range that func is approximated over
 *  \param lowerlimit lower limit of the range that func is approximated over
 *  \param ncoeff number of coefficients to use in the Chebyshev approximation
 *  \param *func pointer to the user-defined function func which is to be approximated
 *
 *  This routine is basically a wrapper around CChebyshevApprox1D::Setup.
 */

/** 
 *  \fn CChebyshevApprox1D::~CChebyshevApprox1D(void)
 *  \brief Plain vanilla destructor.  STL vector class deallocates 
 *         the coefficients for me.
 */

/** 
 *  \fn CChebyshevApprox1D::SetParameters(double, double, int)
 *  \brief Little function to set the upper & lower limits and 
 *        number of coefficients
 *  \param upperlimit upper limit of the range that func is approximated over
 *  \param lowerlimit lower limit of the range that func is approximated over
 *  \param ncoeff number of coefficients to use in the Chebyshev approximation
 */
void CChebyshevApprox1D::SetParameters(double upperlimit,double lowerlimit,  
   int ncoeff) 
{
   // Set limits
   uplim=upperlimit;
   lolim=lowerlimit;
   // Set number of coeffs.
   nc=ncoeff;
}

void CChebyshevApprox1D::Setup(double upperlimit, double lowerlimit, int ncoeff, 
   Array1D<double> coeffs)
{
   // Reset parameters
   this->SetParameters(upperlimit, lowerlimit, ncoeff); 
   // Reallocate memory for coeff. array if needed
   if (c.dim() != nc) {
      Array1D<double> _ctmp(nc);
      _ctmp=0.0;
      c=_ctmp;
   }   
   // Set coeff array
   if ((coeffs.dim() != nc)&&(!quiet)) cerr << "wrong number of coefficients!\n";
   c=coeffs;
};

/** 
 *  \fn CChebyshevApprox1D::Setup(double upperlimit, double lowerlimit, int ncoeff, double (*func)(double, void*), void *pars)
 *  \brief Routine to setup the 1d Chebyshev approximation.
 *  \param upperlimit upper limit of the range that func is approximated over
 *  \param lowerlimit lower limit of the range that func is approximated over
 *  \param ncoeff number of coefficients to use in the Chebyshev approximation
 *  \param *func pointer to the user-defined function func which is to be approximated
 *
 *  This function computes the coefficients of the Chebyshev approximation to 
 *  \c func \c using the well known procedure explained in Numerical Recipes
 *  pp. 191-193.  In detail, the coefficients are computed by:
 *  \f[
 *     c_j = \frac{2}{\mbox{nc}}\sum_{k=0}^{\mbox{nc}-1}
 *     \mbox{func}\left[\cos{\left(\frac{\pi (k+\frac{1}{2})}{\mbox{nc}}\right)}\right] 
 *     \cos\left(\frac{\pi j (k+\frac{1}{2})}{\mbox{nc}}\right) 
 *  \f]
 */
void CChebyshevApprox1D::Setup(double upperlimit, double lowerlimit, 
   int ncoeff, double (*func)(double, void*), void *pars)
{
   //Reset parameters
   this->SetParameters(upperlimit,lowerlimit,ncoeff);
   // Reallocate memory for coeff. array if needed
   if (c.dim() != nc) {
      Array1D<double> _ctmp(nc);
      _ctmp=0.0;
      c=_ctmp;
   }   
   // Set the center and size of the interval we're approximating over --
   // we'll use these to rescale the problem to the interval (0,1)
   double center = 0.5*(uplim+lolim);
   double halfrange = 0.5*(uplim-lolim);
   // Compute the value of the function we're approximating at the 
   // magic points required by Num. Rec. (see Eq. 5.8.7)
   Array1D<double> ctmp(nc);
   for (int i=0;i<nc;i++){
      double y=this->ColocPoint(i);  
      ctmp[i]=func(y*halfrange+center, pars);
   }
   // Now compute the coefficients themselves
   double factor=2.0/nc;
   for(int i=0;i<nc;i++){
      double sum=0.0;
      for (int j=0;j<nc;j++)
         sum += ctmp[j]*cos(PI*i*(j+0.5)/nc);
      c[i] = factor*sum;
   }
}

/**
 *   \fn void CChebyshevApprox1D::SetupFit(double upperlimit, double lowerlimit, int ncoeff, const Vector<double> &controlpts, const Vector<double> &ftnvals)
 *   \brief Routine to fit a Chebyshev polynomial to a table of data.
 *   \param upperlimit upper limit of the fitting region.
 *   \param lowerlimit lower limit of the fitting region.
 *   \param ncoeff number of coefficients to use in the approximation of \f$ f \f$
 *   \param controlpts function arguments \f$ x_i \f$
 *   \param ftnvals value of the function at \f$ x_i \f$, \f$ f_i=f(x_i) \f$
 *
 *    This routine attempts to solve the linear equation 
 *    \f[
 *       f_i = \sum_{j=0}^{\mbox{\tt ncoeff}} c_j T_j(x_i) -\frac{1}{2} c_0
 *    \f]
 *    for the coefficients of the Chebyshev approximation of \f$ f, \f$ c .
 *    It does this by computing the matrix of Chebyshev polynomials evaluated 
 *    at the control points \f$ T_j(x_i) \f$ and then using the Singular Value
 *    Decomposition to make an "inverse" of this matrix.  With this inverse,
 *    the routine solves for the coefficients.  This routine uses
 *    the various SVD routines from jama_svd.h and tnt_linalg.h.
 * 
 *    Note that in the guts of this routine, x and y get rescaled so that they
 *    are between (-1,1).
 */
void CChebyshevApprox1D::SetupFit(double upperlimit, double lowerlimit, 
   int ncoeff, const Array1D<double> &orig_controlpts, 
   const Array1D<double> &orig_ftnvals)
{
   // Set up working vectors
   Array1D<double> controlpts,ftnvals;
   controlpts=orig_controlpts;
   ftnvals=orig_ftnvals;

   // Check if control points and ftn. vals. have same dim
   if ((controlpts.dim() != ftnvals.dim())&&(!quiet)) {
      cerr << "Dimensionality of control point and function value" << 
         "Vectors are different!\n";
   }
   // Check if control points are in the interval
   if ((controlpts[0]<lowerlimit)||(controlpts[controlpts.dim()-1]>upperlimit))
   {
      if (!quiet) {
        cerr << "Control points don't fit in your fitting interval!!\n";
        cerr << "Dammit, i have to do some data pruning\n";
      }
      // count number of control points in interval
      int numgoodpts=0;
      for (int i=0;i<controlpts.dim();i++){
        if ((controlpts[i]>=lowerlimit)&&(controlpts[i]<=upperlimit)) 
          numgoodpts+=1;
      }
      if (!quiet) {cerr << "Using only "<<numgoodpts<<" points in fit"<<endl;}
      // make temp vectors with correct dimensions
      Array1D<double> controltmp(numgoodpts);
      Array1D<double> ftntmp(numgoodpts);
      // copy only the good points
      int igood=0;
      for (int i=0;i<controlpts.dim();i++){
        if ((controlpts[i]>=lowerlimit)&&(controlpts[i]<=upperlimit)) {
            controltmp[igood]=controlpts[i];
            ftntmp[igood]=ftnvals[i];
            igood+=1;
          }
      }
      // replace the real control pt. vector and ftn. val. vector w/ temps
      controlpts=controltmp;
      ftnvals=ftntmp;
   }
   // Check if there's enough data to constrain the approx.
   if ((controlpts.dim() < ncoeff)&&(!quiet)) {
      cerr << "You don't have enough control points to constrain the "<<
         "requested number of coefficients!\n";
   }

   // Setup matrix we will invert
   int M=controlpts.dim(), N=ncoeff;
   Array2D<double> T(M,N);  

   // Reset parameters
   this->SetParameters(upperlimit,lowerlimit,ncoeff);

   // Reallocate memory for coeff. array if needed (nc should be set by now)
   if (c.dim() != nc) {
      Array1D<double> _ctmp(nc);
      _ctmp=0.0;
      c=_ctmp;
   }   

   // Compute center and range of the fit region
   double midpt=0.5*(upperlimit+lowerlimit);
   double halfrange=0.5*(upperlimit-lowerlimit);
   
   // Compute the matrix we want to invert
   {
      double x;int i,j;
      for (i=0;i<M;i++){
         x=(controlpts[i]-midpt)/halfrange;
         T[i][0]=ChebyshevPolynomial(0,x)-0.5;
         for (j=1;j<N;j++){
            T[i][j]=ChebyshevPolynomial(j,x);
         }
      }
   }
   
   // Solve for the coefficient matrix using SVD decomposition
      c=svdinvert(T)*ftnvals;
   
   return;
};

/**
 *  \fn double CChebyshevApprox1D::ColocPoint(int i) const
 *  \brief Routine to return the Colocations points of the Chebyshev
 *         approximation.
 *  \param i index of the Colocation point you are after
 *   
 *  The Colocation points of the Chebyshev approximation are given by    
 *  \f[
 *     x_k=\cos{\left(\pi(k+\frac{1}{2})/\mbox{nc}\right)},   
 *     k=0,1,...,\mbox{nc}
 *  \f]
 */
double CChebyshevApprox1D::ColocPoint(int i) const{
   return cos(PI*(i+0.5)/nc);
}

/** 
 *  \fn double CChebyshevApprox1D::Val(double x, int m) const
 *  \brief Routine to return the approximated value of the 1d Chebyshev 
 *         approximation using the first m coefficients
 *  \param x arguement of func that you which to evaluate the approximation at
 *  \param m number of coefficients you wish to use in the approximation
 *
 *  The Chebyshev approximation to the function \c func \c is given by 
 *  \f[ 
 *     \mbox{func}(x)=\left[\sum_{k=0}^{m-1} c_k T_k(y) \right]-\frac{1}{2}c_0
 *  \f]
 *  where \f$ T(y) \f$ is the Chebyshev polynomial, \f$ c_k=\f$ \c c [k] \c is the 
 *  coefficient of the expansion and \f$ y\f$ is the value of \f$ x\f$ rescaled 
 *  by the upper and lower limits of the approximation:
 *  \f[ 
 *     y=\frac{x-\frac{1}{2}(\mbox{uplim}+\mbox{lolim})}
 *        {\frac{1}{2}(\mbox{uplim}-\mbox{lolim})}
 *  \f] 
 *  The routine doesn't actually do the sum this way -- it actually uses the 
 *  Clenshaw recurrance relation.  See Numerical Recipes pp. 191-193 for more 
 *  details.
 */
double CChebyshevApprox1D::Val(double x, int m) const
{
   // Check the validity of the arguments
   if ((x>uplim)||(x<lolim)) {
      if (!quiet){
        cerr << "x="<<x<<" is out of range in CChebyshevApprox1D"<<endl;
        cerr << "lolim="<<lolim<<"\n";
        cerr << "uplim="<<uplim<<"\n";
      }
      return 0.0;
   }
   if ((m<=0)||(m>nc)) {
      if (!quiet){
        cerr << "m="<<m<<" is an illegal coeff. index in CChebyshevApprox1D"<<endl;
        cerr << "nc="<<nc<<"\n";
      }
      return 0.0;
   }
   // Initialize variables   
   double y,y2,d=0.0,dd=0.0,sv;
   // Rescale the x that we pass in to live on the interval (-1,1)
   y2=2.0*(y=(2.0*x-lolim-uplim)/(uplim-lolim));  
   // Compute the actual approximate value of the function
   for (int j=m-1;j>=1;j--){
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
   }   
   return y*d-dd+0.5*c[0];
}

/**
 *   \fn CChebyshevApprox1D CChebyshevApprox1D::Derivative(void) const
 */
CChebyshevApprox1D CChebyshevApprox1D::Derivative(void) const
{
	int j;
	double con;
   CChebyshevApprox1D deriv;  
   deriv.SetParameters(uplim,lolim,nc); //set up the using the original function
   {
      Array1D<double> _ctmp(deriv.nc);
      _ctmp=0.0;
      deriv.c=_ctmp;
   }

   /* perform the derivative using Numerical Recipes trick on p.195 */
	deriv.c[nc-1]=0.0;
	deriv.c[nc-2]=2*(nc-1)*c[nc-1];
	for (j=nc-3;j>=0;j--) deriv.c[j]=deriv.c[j+2]+2*(j+1)*c[j+1];
	con=2.0/(uplim-lolim);
	for (j=0;j<nc;j++) deriv.c[j] *= con;

   return deriv;
};

/**
 *   \fn CChebyshevApprox1D CChebyshevApprox1D::Integral(void) const
 */
CChebyshevApprox1D CChebyshevApprox1D::Integral(void) const
{
	int j;
	double sum=0.0,fac=1.0,con;
   CChebyshevApprox1D integral;  
   integral.SetParameters(uplim,lolim,nc); //set up the using the original function
   {
      Array1D<double> _ctmp(integral.nc);
      _ctmp=0.0;
      integral.c=_ctmp;
   }

   /* perform the inetgral using Numerical Recipes trick on p.195 */
	con=0.25*(uplim-lolim);
	for (j=1;j<=nc-2;j++) {
		integral.c[j]=con*(c[j-1]-c[j+1])/j;
		sum += fac*integral.c[j];
		fac = -fac;
	}
	integral.c[nc-1]=con*c[nc-2]/(nc-1);
	sum += fac*integral.c[nc-1];
	integral.c[0]=2.0*sum;

   /* now we must resize the derivative object */
   return integral;
};
 
/**
 *  \fn int CChebyshevApprox1D::Ncoeffs(void) const
 *  \brief Returns the number of coefficients in the approximation.
 */

/**
 *  \fn double CChebyshevApprox1D::Upperlimit(void) const
 *  \brief Returns the upper limit of the range of \c func \c to be approximated.
 */

/**
 *  \fn double CChebyshevApprox1D::Lowerlimit(void) const
 *  \brief Returns the lower limit of the range of \c func \c to be approximated.
 */

/**
 *  \fn double CChebyshevApprox1D::Coeff(int i) const
 *  \brief Returns the \c i \c th coeffient of the approximation.
 *  \param i the index of the coefficient
 */

/**
 *  \fn Vector<double> CChebyshevApprox1D::CoeffVector<double>(void) const
 *  \brief Returns the entire \c c \c coeffient vector of the approximation.
 */

/**
 *  \fn double CChebyshevApprox1D::Val(double x) const
 *  \brief Returns the value of the approximation at x, using all \c nc \c
 *         coefficients.
 *  \param x the value to evaluate the approximation at.
 *
 *  This routine is a wrapper around CChebyshevApprox1D::Val.
 */

/**
 *  \fn double CChebyshevApprox1D::operator()(double x, int m) const
 *  \brief Returns the value of the approximation at x, the first \c m \c
 *         coefficients.
 *  \param x the value to evaluate the approximation at.
 *  \param m the number of coefficients to use in the approximation.
 *
 *  This routine is a wrapper around CChebyshevApprox1D::Val.
 */

/**
 *  \fn double CChebyshevApprox1D::operator()(double x) const
 *  \brief Returns the value of the approximation at x, using all \c nc \c
 *         coefficients.
 *  \param x the value to evaluate the approximation at.
 *
 *  This routine is a wrapper around CChebyshevApprox1D::Val.
 */

/**
 *   \fn CChebyshevApprox1D CChebyshevApprox1D::Copy(const CChebyshevApprox1D &A);
 */
CChebyshevApprox1D CChebyshevApprox1D::Copy(const CChebyshevApprox1D &A){
   if (&A != this) {
      uplim=A.uplim;
      lolim=A.lolim;
      quiet=A.quiet;
      nc=A.nc;
      c=A.c;    // let Vector do the dirty work of reallocating & resizing
   }
   return *this;
};

/**
 *  \fn ostream& CChebyshevApprox1D::operator<<(ostream &s, const CChebyshevApprox1D &A)
 */
ostream& operator<<(ostream &s, const CChebyshevApprox1D &A)
{
   s<<A.uplim<<" "<<A.lolim<<endl<<A.c;
   return s;
}

/**
 *  \fn istream& CChebyshevApprox1D::operator>>(istream &s, CChebyshevApprox1D &A)
 */
istream& operator>>(istream &s, CChebyshevApprox1D &A)
{
   s>>A.uplim>>A.lolim>>A.c;
   return s;
}

/**
 *  \fn Array1D< double > CChebyshevApprox1D::c<double>
 *
 *  Vector of coefficients of the Chebyshev polynomial approximation
 */

/**
 *  \fn CChebyshevApprox1D::lolim 
 *
 *  Lower limit of the fitting interval of the Chebyshev polynomial 
 *  approximation
 */

/**
 *  \fn double CChebyshevApprox1D::uplim
 *  
 *  Upper limit of the fitting interval of the Chebyshev polynomial 
 *  approximation
 */

/**
 *  \fn nt CChebyshevApprox1D::nc
 *  
 *  Number of coefficients in the Chebyshev polynomial approximation
 */

//----------------------------------------------------------//
//                                                          //
// Main class definition file for CChebyshevApprox2D.       //
// Documentation is written using Doxygen.                  //
//                                                          //
//----------------------------------------------------------//

/** 
 *  \fn CChebyshevApprox2D::CChebyshevApprox2D(void) 
 *  \brief Plain vanilla constructor for a 2D approximation
 */
CChebyshevApprox2D::CChebyshevApprox2D(void) 
{
   nc_x = 0;
   nc_y = 0;
   uplim_x = 0.0;
   lolim_x = 0.0;
   uplim_y = 0.0;
   lolim_y = 0.0;
}

/** 
 *  \fn CChebyshevApprox2D::~CChebyshevApprox2D(void) 
 *  \brief Plain vanilla destructor for a 2D approximation
 */

/** 
 *  \fn CChebyshevApprox2D::CChebyshevApprox2D(double upperlimit_x, double lowerlimit_x, double upperlimit_y, double lowerlimit_y, int ncoeff_x, int ncoeff_y, double(* func, void*)(double, double),void * pars) 
 *  \brief Basically a wrapper for CChebyshevApprox2D::Setup 
 */

/** 
 *  \fn CChebyshevApprox2D CChebyshevApprox2D::Copy(const CChebyshevApprox2D &A) 
 */
CChebyshevApprox2D CChebyshevApprox2D::Copy(const CChebyshevApprox2D &A)
{
   if (&A != this) {
      uplim_x=A.uplim_x;
      lolim_x=A.lolim_x;
      uplim_y=A.uplim_y;
      lolim_y=A.lolim_y;
      nc_x=A.nc_x;
      nc_y=A.nc_y;
      c=A.c;    // let Matrix do the dirty work of reallocating & resizing
   }
   return *this;
}

/** 
 *  \fn void CChebyshevApprox2D::SetParameters(double upperlimit_x, double lowerlimit_x, double upperlimit_y, double lowerlimit_y, int ncoeff_x, int ncoeff_y) 
 *  \brief Initializes some parameters of 2D approximation
 *  \param upperlimit_x upper limit of the interval in the x direction
 *  \param lowerlimit_x lower limit of the interval in the x direction
 *  \param upperlimit_y upper limit of the interval in the y direction
 *  \param lowerlimit_y lower limit of the interval in the y direction
 *  \param ncoeff_x number of coefficients to use in the x direction
 *  \param ncoeff_y number of coefficients to use in the y direction
 */
void CChebyshevApprox2D::SetParameters(
         double upperlimit_x, double lowerlimit_x, 
         double upperlimit_y, double lowerlimit_y, 
         int ncoeff_x, const int ncoeff_y) 
{
   nc_x = ncoeff_x;
   nc_y = ncoeff_y;
   uplim_x = upperlimit_x;
   lolim_x = lowerlimit_x;
   uplim_y = upperlimit_y;
   lolim_y = lowerlimit_y;
}

/**
 *  \fn double CChebyshevApprox2D::ColocPoint_x(int i) const
 *  \brief Routine to return the x Colocations points of the 2d Chebyshev
 *         approximation.
 *  \param i index of the Colocation point you are after
 *   
 *  The Colocation points of the Chebyshev approximation are given by    
 *  \f[
 *     x_k=\cos{\left(\pi(k+\frac{1}{2})/\mbox{\tt nc\_x}\right)},   
 *     k=0,1,...,\mbox{\tt nc\_x}
 *  \f]
 */
double CChebyshevApprox2D::ColocPoint_x(int i) const{
   return cos(PI*(i+0.5)/nc_x);
}

/**
 *  \fn double CChebyshevApprox2D::ColocPoint_y(int i) const
 *  \brief Routine to return the y Colocations points of the 2d Chebyshev
 *         approximation.
 *  \param i index of the Colocation point you are after
 *   
 *  The Colocation points of the Chebyshev approximation are given by    
 *  \f[
 *     y_k=\cos{\left(\pi(k+\frac{1}{2})/\mbox{\tt nc\_y}\right)},   
 *     k=0,1,...,\mbox{\tt nc\_y}
 *  \f]
 */
double CChebyshevApprox2D::ColocPoint_y(int i) const{
   return cos(PI*(i+0.5)/nc_y);
}

/** 
 *  \fn void CChebyshevApprox2D::Setup(double upperlimit_x, double lowerlimit_x, double upperlimit_y, double lowerlimit_y, int ncoeff_x, int ncoeff_y, double (*func)(double,double, void* ),void *pars)
 *  \brief Routine to setup the 2d Chebyshev approximation
 *  \param upperlimit_x upper limit of approximation in x direction
 *  \param lowerlimit_x lower limit of approximation in x direction
 *  \param upperlimit_y upper limit of approximation in y direction
 *  \param lowerlimit_y lower limit of approximation in y direction
 *  \param ncoeff_x number of coefficients in x direction
 *  \param ncoeff_y number of coefficients in y direction
 *  \param *func pointer to function to be approximated
 *
 *  This routine sets up the 2d Chebyshev approximation by computing the 
 *  coefficient matrix \c c \c using the 2d generalization of the scheme 
 *  used in the 1d approximation:
 *  \f[
 *     c_{ij}=\frac{4}{\mbox{nc\_x}*\mbox{nx\_y}}\sum_{k=0}^{\mbox{nc\_x}}
 *     \sum_{l=0}^{\mbox{nc\_y}} 
 *     f\left(\cos{\left(\frac{\pi(k+\frac{1}{2})}{\mbox{\tt nc\_x}}\right)},
 *            \cos{\left(\frac{\pi(l+\frac{1}{2})}{\mbox{\tt nc\_y}}\right)}
 *      \right)
 *     \cos{\left(\frac{\pi i(k+\frac{1}{2})}{\mbox{\tt nc\_x}}\right)}
 *     \cos{\left(\frac{\pi j(l+\frac{1}{2})}{\mbox{\tt nc\_y}}\right)}
 *  \f]
 */
void CChebyshevApprox2D::Setup(
   double upperlimit_x, double lowerlimit_x, 
   double upperlimit_y, double lowerlimit_y,
   int ncoeff_x, int ncoeff_y, 
   double (*func)(double,double,void*),void* pars) 
{
   //clog << "Computing coefficients for Chebyshev Approx.\n";
   //Reset parameters
   (*this).SetParameters(
      upperlimit_x,lowerlimit_x,
      upperlimit_y,lowerlimit_y,
      ncoeff_x,ncoeff_y);
   // Reallocate memory for coeff. matrix
   if ((c.dim1() != nc_x)||(c.dim2() != nc_y)) {
      Array2D<double> _ctmp(nc_x,nc_y);
      _ctmp=0.0;
      c=_ctmp;
   }   
   // Set the center and size of the interval we're approximating over --
   // we'll use these to rescale the problem to the interval (0,1)
   double center_x = 0.5*(uplim_x+lolim_x);
   double halfrange_x = 0.5*(uplim_x-lolim_x);
   double center_y = 0.5*(uplim_y+lolim_y);
   double halfrange_y = 0.5*(uplim_y-lolim_y);
   // Compute the value of the function at the (rescaled) Colocation points
   Array2D<double> ctmp(nc_x,nc_y); 
   for (int i=0;i<nc_x;i++){
      for (int j=0;j<nc_y;j++){
         ctmp[i][j]=func(this->ColocPoint_x(i)*halfrange_x+center_x,
                         this->ColocPoint_y(j)*halfrange_y+center_y,pars);
      }
   }
   // Compute the coefficients
   double factor=4.0/nc_x/nc_y;
   for (int i=0;i<nc_x;i++){
      for (int j=0;j<nc_y;j++){
         double sum=0.0;
         for (int k=0;k<nc_x;k++){
            for (int l=0;l<nc_y;l++){
               sum+=ctmp[k][l]*cos(PI*i*(k+0.5)/nc_x)*cos(PI*j*(l+0.5)/nc_y);
            }
         }
         c[i][j]=factor*sum;
      }
   }

}

/**
 *  \fn void CChebyshevApprox2D::SetupFit(double upperlimit_x, double lowerlimit_x, double upperlimit_y, double lowerlimit_y, int ncoeff_x, int ncoeff_y, const Vector<double> &x_controlpts, const Vector<double> &y_controlpts, const Matrix<double> &ftnvals)
 *  \brief This routine sets up the 2D Chebyshev polynomial fit of a tabulated function. 
 *  \param upperlimit_x upper limit of approximation in x direction
 *  \param lowerlimit_x lower limit of approximation in x direction
 *  \param upperlimit_y upper limit of approximation in y direction
 *  \param lowerlimit_y lower limit of approximation in y direction
 *  \param ncoeff_x number of coefficients in x direction
 *  \param ncoeff_y number of coefficients in y direction
 *  \param x_controlpts Vector of control points to be used in the x direction, \f$ x_i. \f$
 *  \param y_controlpts Vector of control points to be used in the y direction, \f$ y_j. \f$
 *  \param ftnvals Array2D of function values to fit, \f$ f_{ij}=f(x_i,y_j). \f$
 *
 *    This routine attempts to solve the linear equation 
 *    \f[
 *       f_{ij} = f(x_i,y_j) = \sum_{k=0}^{\mbox{\tt ncoeff\_x}}  
 *                             \sum_{\ell=0}^{\mbox{\tt ncoeff\_y}}
 *                             c_{k\ell} 
 *                             \left(T_k(\tilde{x}_i) -\frac{1}{2} \delta_{0k}\right)
 *                             \left(T_\ell(\tilde{y}_j) -\frac{1}{2} \delta_{0\ell}\right)
 *    \f]
 *    for the coefficients of the Chebyshev approximation of \f$ f, \f$ c .
 *    Note that the tilde means that we rescale the arguments in order that the
 *    range goes from (-1,1).
 *
 *    This routine works by computing the matrix of Chebyshev polynomials evaluated 
 *    at the control points (e.g. \f$ T_j(x_i) \f$ ) and then using the Singular Value
 *    Decomposition to make an "inverse" of this matrix.  With this inverse,
 *    the routine solves for the coefficients.  This routine uses
 *    the various SVD routines from tnt_linalg.h and jama_svd.h.
 *
 */
void CChebyshevApprox2D::SetupFit(double upperlimit_x, double lowerlimit_x, 
   double upperlimit_y, double lowerlimit_y, int ncoeff_x, int ncoeff_y,
   const Array1D<double> &x_controlpts, const Array1D<double> &y_controlpts, 
   const Array2D<double> &ftnvals)
{
   //clog << "Computing coefficients for Chebyshev Approx.\n";
   // Check if everything dimensioned OK
   if (x_controlpts.dim() != ftnvals.dim1()) cerr << 
      "Dimensionality of x control point and function value " << 
      "Vectors are different!\n";
   if (y_controlpts.dim() != ftnvals.dim2()) cerr << 
      "Dimensionality of y control point and function value " << 
      "Vectors are different!\n";
   if (x_controlpts.dim() < ncoeff_x) cerr << 
      "You don't have enough x control points to constrain the "<<
      "requested number of coefficients!\n";
   if (y_controlpts.dim() < ncoeff_y) cerr << 
      "You don't have enough y control points to constrain the "<<
      "requested number of coefficients!\n";
   // Check if control points are in the interval
   if ((x_controlpts[0]<lowerlimit_x)||
      (x_controlpts[x_controlpts.dim()-1]>upperlimit_x)) 
      cerr << "x Control points don't fit in your fitting interval!!\n";
   if ((y_controlpts[0]<lowerlimit_y)||
      (y_controlpts[y_controlpts.dim()-1]>upperlimit_y)) 
      cerr << "y Control points don't fit in your fitting interval!!\n";

   //Reset parameters
   (*this).SetParameters(
      upperlimit_x,lowerlimit_x,
      upperlimit_y,lowerlimit_y,
      ncoeff_x,ncoeff_y);
   // Reallocate memory for coeff. matrix
   if ((c.dim1() != nc_x)||(c.dim2() != nc_y)) {
      Array2D<double> _ctmp(nc_x,nc_y);
      _ctmp=0.0;
      c=_ctmp;
   }   
   
   // Set the center and size of the interval we're approximating over --
   // we'll use these to rescale the problem to the interval (0,1)
   double center_x = 0.5*(uplim_x+lolim_x);
   double halfrange_x = 0.5*(uplim_x-lolim_x);
   double center_y = 0.5*(uplim_y+lolim_y);
   double halfrange_y = 0.5*(uplim_y-lolim_y);
   

   clog << "   Building fit matrix\n";
   // Compute the matrix we want to invert
   int Mpts=x_controlpts.dim()*y_controlpts.dim(), Ncoeffs=ncoeff_x*ncoeff_y;  
   Array2D<double> T(Mpts,Ncoeffs); 
   {
      double x,y;
      double xcoeff,ycoeff;
      int ipt=0;
      for (int i=0;i<x_controlpts.dim();i++){
         x=(x_controlpts[i]-center_x)/halfrange_x;
         for (int j=0;j<y_controlpts.dim();j++){
            y=(y_controlpts[j]-center_y)/halfrange_y;
            int icoeff=0;
            for (int k=0;k<ncoeff_x;k++){
               if ( k == 0 ) {xcoeff=ChebyshevPolynomial(0,x)-0.5;} 
               else {xcoeff=ChebyshevPolynomial(k,x);} 
               for (int l=0;l<ncoeff_y;l++){
                  if ( l == 0 ) {ycoeff=ChebyshevPolynomial(0,y)-0.5;} 
                  else {ycoeff=ChebyshevPolynomial(l,y);} 
                  T[ipt][icoeff]=xcoeff*ycoeff;
                  icoeff+=1;
               }
            }
            ipt+=1;
         }
      }
   }

   clog << "   Inverting to find Matrix of coefficients\n";
   // Pack ftnvals into ftnvalstmp
   Array1D<double> ftnvalstmp(Mpts);
   {
      int ipt=0;
      for (int i=0;i<x_controlpts.dim();i++){
         for (int j=0;j<y_controlpts.dim();j++){
            ftnvalstmp[ipt]=ftnvals[i][j];
            ipt+=1;
         }
      }
   }
   
   // Solve for the coefficient matrix using SVD decomposition
   Array1D<double> ctmp(Ncoeffs);   
   ctmp=svdinvert(T)*ftnvalstmp;   

   // Unpack ctmp into c
   {
      int icoeff=0;
      for (int k=0;k<ncoeff_x;k++){
         for (int l=0;l<ncoeff_y;l++){
            c[k][l]=ctmp[icoeff];
            icoeff+=1;
         }
      }
   }
   clog << "Coefficients for Chebyshev Approx. computed!\n";
}

/**
 *  \fn double CChebyshevApprox2D::Val(double x, double y, int m_x, int m_y) const
 *  \brief Routine to return the approximated value of the Chebyshev 
 *         approximation using the first m_x coefficients in the x direction 
 *         and the first m_y coefficients in the y direction
 *  \param x  x coordinate you want the approximation evaluated at
 *  \param y  y coordinate you want the approximation evaluated at
 *  \param m_x number of coefficients to keep in the x direction
 *  \param m_y number of coefficients to keep in the y direction
 *
 *  This function returns the approximated value of the Chebyshev using
 *  a generalization of the procedure advocated in Numerical Recipes 
 *  pp. 191-193.  Basically, we write the approximation of \c func(x,y) \c in 
 *  the x direction as 
 *  \f[ 
 *     \mbox{func}(x,y)=\left[\sum_{k=0}^{m_x-1} c_k(y) T_k(\tilde{x})
 *     \right]-\frac{1}{2}c_0(y)
 *  \f]
 *  where \f$ T(\tilde{x}) \f$ is the Chebyshev polynomial, \f$ c_k(y)\f$ is the 
 *  coefficient of the expansion as a function of \f$ y \f$ and \f$ \tilde{x}\f$ 
 *  is the value of \f$ x\f$ rescaled 
 *  by the upper and lower limits of the approximation:
 *  \f[ 
 *     \tilde{x}=\frac{x-\frac{1}{2}(\mbox{uplim\_x}+\mbox{lolim\_x})}
 *        {\frac{1}{2}(\mbox{uplim\_x}-\mbox{lolim\_x})}
 *  \f] 
 *  The routine doesn't actually do the sum this way -- it actually uses the 
 *  Clenshaw recurrance relation.  See Numerical Recipes pp. 191-193 for more 
 *  details.
 *
 *  How do we compute \f$ c_k(y)\f$?  Well, basically the same way...  we make
 *  a Chebyshev approximation of it:
 *  \f[
 *     c_k(y)=\left[\sum_{l=0}^{m_y-1} c_{kl} T_l(\tilde{y})
 *     \right]-\frac{1}{2}c_{k0}
 *  \f]
 *  Here the terms mean the same thing as above, but noting that the 
 *  coefficient matrix \f$ c_{kl}=\f$ \c c [k][l], \c which we compute 
 *  above in CChebyshevApprox2D::Setup.
 */
double CChebyshevApprox2D::Val(double x, double y, int m_x, int m_y) const
{
   // Check the validity of the function arguments
   if ((x > uplim_x) || (x < lolim_x)) cout << "x argument outside of region" <<
      " of validity of Chebyshev approximation. " << endl;
   if ((y > uplim_y) || (y < lolim_y)) cout << "y argument outside of region" <<
      " of validity of Chebyshev approximation. " << endl;
   if ((m_x > nc_x) || (m_x < 1)) cout << "illegal m_x requested" << endl;
   if ((m_y > nc_y) || (m_y < 1)) cout << "illegal m_y requested" << endl;

   // Initialize variables   
   double xt2,xt,yt2,yt,d,dd,sv;
   // Get the arguements, but rescaled to -1,1 for the approx.
   xt2=2.0*(xt=(2.0*x-lolim_x-uplim_x)/(uplim_x-lolim_x));
   yt2=2.0*(yt=(2.0*y-lolim_y-uplim_y)/(uplim_y-lolim_y));
   // The coefficients of the chebyshev polynomial approximation to the 
   // chebyshev polynomial coeffs.
   Array1D<double> ctmp(nc_x);
   // Compute the approximate value of the coefficient of the chebyshev
   // approximation of the coefficient in the x direction
   for (int i=0;i<nc_x;i++){
      d=dd=0.0;
      for (int j=m_y-1;j>=1;j--){
         sv=d;
         d=yt2*d-dd+c[i][j];
         dd=sv;
      }   
      ctmp[i]=yt*d-dd+0.5*c[i][0];
   }   
   // Compute the actual approximate value of the function
   d=dd=0.0;
   for (int i=m_x-1;i>=1;i--){
      sv=d;
      d=xt2*d-dd+ctmp[i];
      dd=sv;
   }   
   return xt*d-dd+0.5*ctmp[0];
}

/**
 *  \fn ostream& CChebyshevApprox2D::operator<<(ostream &s, const CChebyshevApprox2D &A)
 */
ostream& operator<<(ostream &s, const CChebyshevApprox2D &A)
{
   s<<A.uplim_x<<" "<<A.lolim_x<<endl;
   s<<A.uplim_y<<" "<<A.lolim_y<<endl;
   s<<A.c;
   return s;
}

/**
 *  \fn istream& CChebyshevApprox2D::operator>>(istream &s, CChebyshevApprox2D &A)
 */
istream& operator>>(istream &s, CChebyshevApprox2D &A)
{
   s>>A.uplim_x>>A.lolim_x;
   s>>A.uplim_y>>A.lolim_y;
   s>>A.c;
   return s;
}

