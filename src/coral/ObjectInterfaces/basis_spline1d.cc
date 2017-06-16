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
#include "basis_spline1d.h"
#include <cmath>

#define BSPLINE_EPSILON 1e-12

using namespace std;
using namespace TNT;

// Read/write to parameter map
bool CBasisSpline1d::Read(const parameterMap& m){
    bool success = CGenericSpline1d::Read(m);
    Array1D<double> new_splco(ndata,0.);
    _splco = new_splco;
    return success;
}

/** 
    \brief Get value of Basis Spline at x
    \warning The resulting spline is not left continuous at the left-most knot
*/
double CBasisSpline1d::getValue(double x) const{
    if ( (x<knots[0]) || (x>knots[knots.dim()-1]) ) return 0.0;
    if (std::abs(x-knots[0])<BSPLINE_EPSILON) x+=BSPLINE_EPSILON;
    if (std::abs(x-knots[knots.dim()-1])<BSPLINE_EPSILON) x-=BSPLINE_EPSILON;
    return get_bvalue(data,x,0);
}

/**
    \brief Gets the value of Basis Spline i at x  
    \param i Basis Spline index
    \param x value to evaluate spline at
    \param jderiv order of derivative to take
*/
double CBasisSpline1d::basisFunction(double x, int i, int jderiv) const{
    // Make sure we have legal index
    if  ((i<0) || (i>=ndata)) return 0.0;
    // Make sure this B-spline is non-zero before wasting time evaluating it
    if ((x < knots[i]) || (x > knots[i+spline_degree+1])) return 0.0;
    // Set up dummy array so get just one spline out
    const_cast<CBasisSpline1d*>(this)->_splco[i]=1.0;
    // take care of continuity at edges
    if (std::abs(x-knots[0])<BSPLINE_EPSILON) x+=BSPLINE_EPSILON;
    if (std::abs(x-knots[spline_degree-1])<BSPLINE_EPSILON) x-=BSPLINE_EPSILON;
    double ans=get_bvalue(_splco,x,jderiv);
    // put dummy array back the way it was
    const_cast<CBasisSpline1d*>(this)->_splco[i]=0.0;
    return ans;
}



/**  
   \brief Calculates value at  x  of  jderiv-th derivative of spline from b-repr.
   the spline is taken to be continuous from the right, EXCEPT at the
   rightmost knot, where it is taken to be continuous from the left.

   From  "A Practical Guide to Splines"  by C. De Boor    
    
 \param x the point at which to evaluate .
 \param bcoef b-coefficient sequence, of length  _nc
 \param jderiv integer giving the order of the derivative to be evaluated assumed  to be zero or positive.  
        
 \warning The restriction  k .le. kmax (=20) is imposed arbitrarily by the dimension statement for  aj, dl, dr  below, but is  nowhere  checked  for.
 
 \retval bvalue the value of the (jderiv)-th derivative of  \f$f\f$  at  \f$x\f$ .
  

     The nontrivial knot interval  \f$(t_i,t_{i+1})\f$  containing  \f$x\f$  
  is located with the aid of  interv . The  \f$k\f$  b-coeffs of  \f$f\f$  
  relevant for this interval are then obtained from  _coeffs (or taken to be 
  zero if not explicitly available) and are then differenced  jderiv  times to
  obtain the b-coeffs of  \f$ d^{jderiv}f \f$  relevant for that interval.
  Precisely, with  \f$j = jderiv\f$, we have from x.(12) of the text that

     \f$ d^jf  =  \sum bcoef(.,j)*b(.,k-j,t) \f$

  where
    \f$ bcoef(.,j)  =  \left\{ 
        \begin{array}{ll}
            bcoef(.) & \mbox{if}\; j = 0 \\
            \frac{\displaystyle bcoef(.,j-1) - bcoef(.-1,j-1) }{\displaystyle (t_{.+k-j} - t_{.})/(k-j)} & j > 0
        \end{array}
        \right.
    \f$
    
  Then, we use repeatedly the fact that
    \f$\sum ( a(.)*b_{.,m,t}(x) )  =  \sum ( a(.,x)*b_{.,m-1,t}(x) )\f$
  with
    \f$a(.,x) 
    =\frac{(x-t(.))*a(.)+(t(.+m-1)-x)*a(.-1)}{(x-t(.))+(t(.+m-1)-x)}\f$
  to write  \f$d^j f(x)\f$  eventually as a linear combination of b-splines
  of order  1 , and the coefficient for  \f$b_{i,1,t}(x)\f$  must then be the
  desired number  \f$d^jf(x)\f$. (see x.(17)-(19) of text).
*/
double CBasisSpline1d::get_bvalue(const Array1D<double>& bcoef, double x, int jderiv) const
{
    int i,ilo,imk,j,jc,jcmin,jcmax,jj,kmj,km1,nmi,jdrvp1; // mflag
    const int kmax = 20;
    int n = ndata;
    int k = spline_degree+1;
    Array1D<double> aj(kmax,0.0),dl(kmax,0.0),dr(kmax,0.0);
    double fkmj;
    double bvalue = 0.;
    if (jderiv >= k) return bvalue;
/*  *** Find  i   s.t.   1 .le. i .lt. n+k   and   t(i) .lt. t(i+1)   and
      t(i) .le. x .lt. t(i+1) . If no such i can be found,  x  lies
      outside the support of  the spline  f , hence  bvalue = 0.
      (The asymmetry in this choice of  i  makes  f  rightcontinuous, except
      at  t(n+k) where it is leftcontinuous.)
*/
    i=getKnotToLeft(x)+1;
    if ( (i < 0)||(i > n+k) ) return bvalue;
    if ( (knots[i-1]>x)||(x>=knots[i]) ) return bvalue;
    
/*  *** if k = 1 (and jderiv = 0), bvalue = _coeffs(i).  */
    km1 = spline_degree; //spline_degree=k-1 by construction
    if (spline_degree < 0) {bvalue = bcoef[i-1]; return bvalue;}

/*
  *** store the k b-spline coefficients relevant for the knot interval
     (t(i),t(i+1)) in aj(1),...,aj(k) and compute dl(j) = x - t(i+1-j),
     dr(j) = t(i+j) - x, j=1,...,k-1 . set any of the aj not obtainable
     from input to zero. set any t.s not obtainable equal to t(1) or
     to t(n+k) appropriately.
*/
    jcmin = 1;
    imk = i - k;
    if (imk >= 0) {
        for(j=1;j<=km1;j++){dl[j-1] = x - knots[i-j];}
    } else {
        jcmin = 1 - imk;
        for(j=1;j<=i;j++){dl[j-1] = x - knots[i-j];}
        for (j=i;j<=km1;j++){
            aj[k-j-1] = 0.;
            dl[j-1] = dl[i-1];
        }
    }
    
    jcmax = k;
    nmi = n - i;
    if (nmi >= 0){
        for (j=1;j<=km1;j++){dr[j-1] = knots[i+j-1] - x;}
    } else {
        jcmax = k + nmi;
        for (j=1;j<=jcmax;j++){dr[j-1] = knots[i+j-1] - x;}
        for (j=jcmax;j<=km1;j++){
            aj[j] = 0.;
            dr[j-1] = dr[jcmax-1];
        }
    }

    for( jc=jcmin;jc<=jcmax;jc++){aj[jc-1] = bcoef[imk + jc-1];}

/*               *** difference the coefficients  jderiv  times.  */
    if (jderiv != 0) {
        for (j=1;j<=jderiv;j++){
            kmj = k-j;
            fkmj = static_cast<double>(kmj);
            ilo = kmj;
            for(jj=1;jj<=kmj;jj++) {
                aj[jj-1] = ((aj[jj] - aj[jj-1])/(dl[ilo-1] + dr[jj-1]))*fkmj;
                ilo = ilo - 1;
            }
        }
    }
/*
  *** compute value at  x  in (t(i),t(i+1)) of jderiv-th derivative,
     given its relevant b-spline coeffs in aj(1),...,aj(k-jderiv).
*/
    if (jderiv != km1) {
        jdrvp1 = jderiv + 1;     
        for (j=jdrvp1;j<=km1;j++){
            kmj = k-j;
            ilo = kmj;
            for (jj=1;jj<=kmj;jj++){
                aj[jj-1]=(aj[jj]*dl[ilo-1]+aj[jj-1]*dr[jj-1])/
                         (dl[ilo-1]+dr[jj-1]);
                ilo = ilo - 1;
            }
        }
    }
    bvalue = aj[0];
    
    return bvalue;
}

//! Tool to redim the knots, data and covmtx and keep dimensions in sync.  Will hose array contents
bool CBasisSpline1d::setDim(int ncoeffs){
    Array1D<double> new_splco(ncoeffs,0.);
    _splco = new_splco;
    return CGenericSpline1d::setDim(ncoeffs);
}
