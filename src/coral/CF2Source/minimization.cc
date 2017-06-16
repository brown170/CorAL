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
#include  <cmath> 
#include  <cstdlib> 
#include  <cstdio> 
#include  "minimization.h" 

using   namespace  std; 

// As an abstract class, this will not be used. Only for demonstration purpose and test. 
CMinimization::CMinimization( int  dimension) 
{ 
	TOL_CG = 1.0e-10; 
	TOL_Brent = 1.0e-12; 
	n = dimension; 
	vec_x = NULL; 
	vec_dx = NULL; 
	linear_trial_x = NULL; 
	linear_dir = NULL; 
	H = NULL; 

	AllocVectors(); 
} 

CMinimization::CMinimization() 
{ 
	//TOL_CG = 1.0e-4; 
	//TOL_Brent = 1.0e-5 ; 
	TOL_CG = 1.0e-6; 
	TOL_Brent = 1.0e-7; 
	n = 1; 
	vec_x = NULL; 
	vec_dx = NULL; 
	linear_trial_x = NULL; 
	linear_dir = NULL; 
	H = NULL; 
} 

bool  CMinimization::AllocVectors() 
{ 
	// Clear previous memory storage 
	if (NULL != vec_x) 
	{ 
		delete []  vec_x; 
	} 
	if (NULL != vec_dx) 
	{ 
		delete []  vec_dx; 
	} 
	if (NULL != linear_trial_x) 
	{ 
		delete []  linear_trial_x; 
	} 
	if (NULL != linear_dir) 
	{ 
		delete []  linear_dir; 
	} 
	if (NULL != H) 
	{ 
		for ( int  i=0;i<n;i++) 
		{ 
			delete []  H[i]; 
		} 
		delete []  H; 
	} 

	// Allocation 
	if (n > 1) 
	{ 
		TagMultiD =  true ; 
		vec_x =  new   double [n];   // x 
		vec_dx =  new   double [n];   // gradient 
		linear_trial_x =  new   double [n];  // trial x 
		linear_dir =  new   double [n];  // search direction 
		H =  new   double  * [n];   // Hessian matrix 
		for ( int  i=0;i<n;i++) 
		{ 
			H[i] =  new   double [n]; 
		} 

	} 
	else   if (1 == n) 
	{ 
		TagMultiD =  false ; 
		vec_x = NULL; 
		vec_dx = NULL; 
		linear_trial_x = NULL; 
		linear_dir = NULL; 
		H = NULL; 
	} 
	else 
	{ 
		printf("Error: please input positive integer value of dimension! input n=%d", n); 
		exit(0); 
	} 

	return  true ; 
} 

bool  CMinimization::SetDimension( int  dimension) 
{ 
	n = dimension; 
	if (AllocVectors() ==  true ) 
	{ 
		return  true ; 
	} 
	return  false ; 
} 

CMinimization::~CMinimization() 
{ 
	if (NULL != vec_x) 
	{ 
		delete []  vec_x; 
	} 
	if (NULL != vec_dx) 
	{ 
		delete []  vec_dx; 
	} 
	if (NULL != linear_trial_x) 
	{ 
		delete []  linear_trial_x; 
	} 
	if (NULL != linear_dir) 
	{ 
		delete []  linear_dir; 
	} 
	if (NULL != H) 
	{ 
		for ( int  i=0;i<n;i++) 
		{ 
			delete []  H[i]; 
		} 
		delete []  H; 
	} 
} 

bool  CMinimization::bracket( double  &ax,  double  &bx,  double  &cx,  double  &fa,  double  &fb,  double  &fc) 
// Given ax and bx, returning three points a, b, and c, where a minimum is between [a,c].If found, return true, otherwise return false. 
{ 
	const  double  GOLD_RATIO = 1.618034;  // Gold ratio in magnifiying intervals 
	const  double  SEARCH_LIMIT = 10.0; //100.0;  // The most far to go in search downhill 

	fa = fn1(ax); 
	fb = fn1(bx); 

	if (fb > fa)  // For convenience, we set the downhill direction to be from a to b. 
	{ 
		swap2(ax, bx); 
		swap2(fa, fb); 
	} 

	cx = bx + GOLD_RATIO * (bx - ax);  // First trial of cx 
	fc = fn1(cx); 

	double  u;  // trial of new point to get minimum bracketed 
	double  fu; 
	double  ulim;  // the maximum search range 
	while (fb > fc) 
	{ 
	//  Inverse parabolic interpolation 
	//u = bx - 0.5 * ((b-a)^2*(fb-fc)-(b-c)^2*(fb-fa))/((b-a)*(fb-fc)-(b-c)*(fb-fa)) 
		u = bx - ((bx - cx) * (bx - cx) * (fb - fa) - (bx - ax) * (bx - ax) * (fb - fc)) / (2.0 * ((bx - cx) * (fb - fa) - (bx - ax) * (fb - fc))); 
		ulim = bx + SEARCH_LIMIT * (cx - bx); 

		if ((bx - u) * (u - cx) > 0.0)  // u is in [b,c] 
		{ 
			fu = fn1(u); 
			if (fu < fc)  // a minimum is in [b,u] 
			{ 
				ax = bx; 
				bx = u; 
				fa = fb; 
				fb = fu; 
				return  true ; 
			} 
			else 
			{ 
				if (fu > fb)  // a minimum is in [a,u] 
				{ 
					cx = u; 
					fc = fu; 
					return  true ; 
				} 
				u = cx + GOLD_RATIO * (cx - bx);  // Parabolic fit fails. Use Gold ratio for u. // fu > fc && fu < fb 
				fu = fn1(u); 
			} 
		} 
		else   if ((cx - u) * (u - ulim) > 0.0)  // u is in [c, ulim] 
		{ 
			fu = fn1(u); 
			if (fu < fc) 
			{ 
				shift3(bx, cx, u, cx + GOLD_RATIO * (cx - bx)); 
				shift3(fb, fc, fu, fn1(u)); 
			} 
		} 
		else   if ((u - ulim) * (ulim - cx) > 0.0)  // u is out of ulim, but we keep it within range. 
		{ 
			u = ulim; 
			fu = fn1(u); 
		} 
		else  // u in [a,b] or (-infy,a]. reject this u, use Gold ratio to set u. 
		{ 
			u = cx + GOLD_RATIO * (cx - bx); 
			fu = fn1(u); 
		} 
		shift3(ax, bx, cx, u); 
		shift3(fa, fb, fc, fu); 
	} 
	// minimum is in [a,c] 
	return  true ; 
} 

double  CMinimization::brent( double  ax,  double  bx,  double  cx,  double  &xmin) 
// One-dimensional search using brent's method 
// a, b, and c are in the order that fb is less than both fa and fc. 
// return minimum value and xmin. 
{ 
	const  int  MAX_ITERATION = 100; 
	const  double  GOLD_RATIO = 2.0 - 1.618034; 

	double  a, b; 

	a = (ax < cx ? ax : cx); 
	b = (ax > cx ? ax : cx);  // reset [a,b] in this order for future convenience. so the input ax, cx can be in any order. 

	double  w, v;  // Temporary storage of previous trial of x 
	double  fw, fv; 
	double  x; 
	double  fx; 
	x = bx;  
	w = bx;  
	v = bx; 
	fx = fn1(x);  
	fw = fx;  
	fv = fx; 
	double  u, fu; 

	double  xm;  // mid-point of [a,b] 

	double  e = 0.0; 
	double  d = 0.0; 
	double  p, q; 
	double  tol1; 

	for ( int  i=0;i<MAX_ITERATION;i++) 
	{ 
		xm = 0.5 * (a + b); 
		tol1 = TOL_Brent * fabs(x) + 1.0e-18 ; 
			if (fabs(x-xm) <= (tol1*2.0-0.5*(b-a)))  // two criteria:1, x->xm; and 2, (b-a)->0 
		{ 
			xmin = x; 
			return  fx;  // find! 
		} 
		if (fabs(e) > tol1) 
		{ 
			q = 2.0 * ((x-v)*(fx-fw) - (x-w)*(fx-fv)); 
			p = (x-v)*(x-v)*(fx-fw) - (x-w)*(x-w)*(fx-fv); 

			if (q > 0.0) p = -p; 
			q = fabs(q);  // reset p, q, so that q>0. but still keep the correct sign before p/q; 
			if (fabs(p) >= fabs(0.5*q*e) 
				|| p <= q*(a-x) 
				|| p >= q*(b-x)) 
			{ 
				e = (x >= xm ? a-x : b-x); 
				d = GOLD_RATIO * e; 
			} 
			else 
			{ 
				e = d; 
				d = p / q;  // parabolic interpolation 
				u = x + d; 
				if  ( u - a < tol1*2.0 || b - u < tol1*2.0) 
				{ 
					d = tol1 / (x-xm) * fabs(x-xm); 
				} 

			} 
		} 
		else 
		{ 
			e = (x >= xm ? a-x : b-x); 
			d = GOLD_RATIO * e; 
		} 

		if (fabs(d) > tol1) 
			u = x + d; 
		else 
			u = x + tol1 / d * fabs(d); 
		fu = fn1(u); 

		if (fu <= fx) 
		{ 
			if (u >= x) 
				a = x; 
			else 
				b = x; 
			shift3(v, w, x, u); 
			shift3(fv, fw, fx, fu); 
		} 
		else 
		{ 
			if (u < x) 
				a = u; 
			else 
				b = u; 
			if (fu <= fw || w == x)  // The book-keeping below is to set up the order of [a,b] = [v, w]. 
			{ 
				v = w;  
				w = u; 
				fv = fw;  
				fw = fu; 
			} 
			else   if (fu <= fv || v == x || v == w) 
			{ 
				v = u; 
				fv = fu; 
			} 
		} 

	} 
	printf("Warning: cann't find the minimum after maximum iteration in brent!\n"); 
	xmin = x; 
	return  fx; 

} 

bool  CMinimization::linear_Min( double  &fmin) 
{ 
	// set up initial searching point. these values should not affect the final result   
	//because, although bracket() guarantees a local minimum, instead of a global minimum,  
	//the conjugate gradient method only requires a local minimum to converge. 
	double  ax = 0.0; 
	double  xx = 1.0e-6; //1.0;  // However when the gradient is very large, the initial guess of xx has to be reset to small enough not to jump out. 

	double  max_dir = 0.0; 
	int  j = 0; 
	for ( int  i=0;i<n;i++) 
	{ 
		if (fabs(linear_dir[i]) > max_dir) 
		{ 
			max_dir = fabs(linear_dir[i]); 
			j = i; 
		} 
	} 
	xx = (max_dir > 1.0) ? 1.0/max_dir : 1.0; 

	// This is very tricky for that, when the gradient of the function at certain points is HUGE(1.0e8), this xx is thus set to be extremely small. 
	// To perform the brent search properly, we need to keep a close eye  at  the criteria of brent's convergence 
	// if(fabs(x-xm) <= (tol1*2.0-0.5*(b-a))) 
	// tol1 = TOL_Brent * fabs(x) + 1.0e-12; here 1.0e-12 is to avoid the disaster when fabs(x) happens to be nearby the machine precision and the convergence is never approachable. 
	// in this case, as x is only 1.0e-9, the brent method converges at a precision of only 1.0e-3(1.0e-12/1.0e-9), which is way below what we might expect. 

	double  bx; 
	double  fa, fb, fx; 
	bracket(ax, xx, bx, fa, fx, fb); 
	printf("bracket done as [%g, %g]\n", ax, bx); 
	double  xmin; 
	fmin = brent(ax, xx, bx, xmin); 
	for ( int  i=0;i<n;i++) 
	{ 
		linear_dir[i] *= xmin;  // calculate the offset of vec_x 
		vec_x[i] += linear_dir[i];  // reset x to the one-dimensional minimum point. 
	} 
	return  true ; 
} 

bool  CMinimization::conjugate_gradient( double  * x,  int  &iteration,  double  &fmin) 
{ 
	for ( int  i=0;i<n;i++)  // copy initial guess to vec_x 
	{ 
		vec_x[i] = x[i];  
	} 
	const  int  MAX_ITERATION = 100; 
	double  * g =  new   double [n]; 
	double  * h =  new   double [n]; 

	double  fx = fn(vec_x); 
	dfn(vec_x);  // the gradient is stored in vec_dx 

	for ( int  i=0;i<n;i++)  // initialize g and h 
	{ 
		g[i] = - vec_dx[i]; 
		h[i] = g[i]; 
		vec_dx[i] = h[i];  
	} 
	double  gg, dgg, gam; 
	for ( int  i=0;i<MAX_ITERATION;i++) 
	{ 
		iteration = i;  // return this when complete minimum search 
		for ( int  ii=0;ii<n;ii++) 
		{ 
			linear_dir[ii] = vec_dx[ii]; 
		} 
		linear_Min(fmin); 

		for ( int  kk=0;kk<n;kk++) 
		{ 
			printf("x[%d] : %g  ", kk, vec_x[kk]); 
		} 
		printf("f(x)=%g after %d iterations\n", fmin, i); 

		if (2.0 * fabs(fmin-fx) <= TOL_CG * (fabs(fmin)+fabs(fx)+1.0e-15))  // Using relative criteria for convergence. 
		// adding 1.0e-15 to deal with the possibility that the function has miminum  value 0.0. 
		{ 
			delete []  g; 
			delete []  h; 
			for ( int  ii=0;ii<n;ii++) 
			{ 
				x[ii] = vec_x[ii]; 
			} 
			return   true ; 
		} 
		fx = fmin;  // reset new starting point of direction search 
		dfn(vec_x);  // recalculate gradient at new point 
		gg = 0.0; 
		dgg = 0.0; 
		for ( int  j=0;j<n;j++) 
		{ 
			gg += g[j]*g[j]; 
	//  dgg += vec_dx[j]*vec_dx[j];  // Fletcher_Reeves. 
			dgg += (vec_dx[j]+g[j])*vec_dx[j];  // Polak-Ribiere. 
		} 
		if (fabs(gg) <= 1.0e-12)  // if the gradient is zero, the minimum is reached. but this won't happen often. 
		{ 
			delete []  g; 
			delete []  h; 
			for ( int  ii=0;ii<n;ii++) 
			{ 
				x[ii] = vec_x[ii]; 
			} 
			return  true ; 
		} 
		gam = dgg / gg;  // the previous direction only used in giving this ratio for conjugate direction. 
		for ( int  j=0;j<n;j++) 
		{ 
			g[j] = - vec_dx[j]; 
			h[j] = g[j] + gam * h[j]; 
			vec_dx[j] = h[j]; 
		} 
	} 
	delete []  g; 
	delete []  h; 
	printf("Cannot converge to minimum after %d iterations with f(x) = %g\n", iteration, fx); 
	return  false ; 
} 

double  CMinimization::fn1( double  x) 
{ 
	if ( false  == TagMultiD) 
	{ 
	// test functions 
		return  sin(0.1 * x); 
	//  return  -x*x*x+100*x*x-100*x; 
	//  return  -0.05*x*x*x+x*x-2*x; 
	} 
	else 
	{ 
	// This is for finding the function's value when variables offset by x along vec_dir direction! 
		for ( int  i=0;i<n;i++) 
		{ 
			linear_trial_x[i] = vec_x[i] + x*linear_dir[i]; 
		} 
		return  fn(linear_trial_x); 
	} 
} 

//bool  CMinimization::dfn(double * x) 
//{ 
//  // test 1 
//  //vec_dx[0] = 2*x[0] + x[1]; 
//  //vec_dx[1] = x[0] + 2*x[1]; 
//  // test 2 
//  //vec_dx[0] = 2*x[0] - 0.5*x[1] + 3; 
//  //vec_dx[1] = -0.5*x[0] + 2*x[1] + 4; 
//  // test 3 
//  //vec_dx[0] = 3*x[0]*x[0]-0.5*x[1]*x[1]+0.5*x[2]*x[2]+0.4*x[0]*x[0]*x[0]; 
//  //vec_dx[1] = -x[0]*x[1]+0.5*x[2]+0.4*x[1]*x[1]*x[1]; 
//  //vec_dx[2] = x[0]*x[2]+0.5*x[1]+0.4*x[2]*x[2]*x[2]; 
//  // test 4 
//  //vec_dx[0] = 3*x[0]*x[0]-0.5*x[1]*x[1]+0.5*x[2]*x[2]+x[3]*x[3]+0.4*x[0]*x[0]*x[0]; 
//  //vec_dx[1] = -x[0]*x[1]+0.5*x[2]+2*x[1]*x[3]+0.4*x[1]*x[1]*x[1]; 
//  //vec_dx[2] = x[0]*x[2]+0.5*x[1]+0.4*x[2]*x[2]*x[2]; 
//  //vec_dx[3] = 3*x[3]*x[3]+2*x[0]*x[3]+x[1]*x[1]+0.4*x[3]*x[3]*x[3]; 
//  return  true; 
//} 

//double  CMinimization::fn(double * x) 
//{ 
//  return  0.0; 
//  // test 1 
//  //return  x[0]*x[0]+x[1]*x[0]+x[1]*x[1]; 
//  // test 2 
//  //return  x[0]*x[0]-0.5*x[1]*x[0]+x[1]*x[1]+3*x[0]+4*x[1]+5; 
//  // test 3 
//  //return  x[0]*x[0]*x[0]-0.5*x[0]*x[1]*x[1]+0.5*x[0]*x[2]*x[2]+0.5*x[1]*x[2]+0.1*(x[2]*x[2]*x[2]*x[2]+x[1]*x[1]*x[1]*x[1]+x[0]*x[0]*x[0]*x[0])+1000; 
//  // test 4 
//  //return  x[0]*x[0]*x[0]-0.5*x[0]*x[1]*x[1]+0.5*x[0]*x[2]*x[2]+0.5*x[1]*x[2]+x[3]*x[3]*x[3] 
//  //+x[0]*x[3]*x[3]+x[1]*x[1]*x[3]+0.1*(x[3]*x[3]*x[3]*x[3]+x[2]*x[2]*x[2]*x[2]+x[1]*x[1]*x[1]*x[1]+x[0]*x[0]*x[0]*x[0])+1000; 
//} 


// for test 
// Before using the below test, one needs to change the pure virtual declarations 
// back into normal function and add the definition of fn and dfn. 
// CMinimization itself is declared as ABSTRACT class. 

//void  main() 
void  dummyfunctionname() 
{ 
	//double  fmin, xmin; 

	//CMinimization  m1(1); 
	//double  ax, bx, cx, fa, fb, fc; 
	//ax = -10; bx = -20; 
	//m1.bracket(ax, bx, cx, fa, fb, fc); 
	//printf("the three points are f(%g) = %g, f(%g) = %g, f(%g) = %g\n", ax, fa, bx, fb, cx, fc); 
	//fmin = m1.brent(ax, bx, cx, xmin); 
	//printf("the minimum is f(%g) = %g\n", xmin, fmin); 

	//CMinimization  m2(2); 
	//double  *  xi = new double[2]; 
	//xi[0] = 1.0; xi[1] = 2.0; 
	//int  iter; 
	//if(m2.conjugate_gradient(xi, iter, fmin)==true) 
	//{ 
	//  printf("find minimum after %d iterations and f(x) = %g\n", iter, fmin); 
	//  printf("x[0] = %g, x[1] = %g\n", xi[0], xi[1]); 
	//  printf("dx[0] = %g, dx[1] = %g\n", m2.vec_dx[0], m2.vec_dx[1]); 
	//} 
	//else 
	//{ 
	//  printf("Conjugate gradient method failed! f(x) = %g\n", fmin); 
	//} 


	//CMinimization  m3(3); 
	//double  *  xi = new double[3]; 
	//xi[0] = -1.0; xi[1] = 2.0e4; xi[2] = -3.0; 
	//int  iter; 
	//if(m3.conjugate_gradient(xi, iter, fmin)==true) 
	//{ 
	//  printf("find minimum after %d iterations and f(x) = %g\n", iter, fmin); 
	//  printf("x[0] = %g, x[1] = %g, x[2] = %g\n", xi[0], xi[1], xi[2]); 
	//  printf("dx[0] = %g, dx[1] = %g, dx[2] = %g\n", m3.vec_dx[0], m3.vec_dx[1], m3.vec_dx[2]); 
	//} 
	//else 
	//{ 
	//  printf("Conjugate gradient method failed! f(x) = %g\n", fmin); 
	//} 
	//delete[]  xi; 

	//CMinimization  m4(4); 
	//double  *  xi = new double[4]; 
	//xi[0] = -1.0; xi[1] = 2.0e4; xi[2] = -3.0; xi[3] = 100; 
	//int  iter; 
	//if(m4.conjugate_gradient(xi, iter, fmin)==true) 
	//{ 
	//  printf("find minimum after %d iterations and f(x) = %g\n", iter, fmin); 
	//  printf("x[0] = %g, x[1] = %g, x[2] = %g, x[3] = %g\n", xi[0], xi[1], xi[2], xi[3]); 
	//  printf("dx[0] = %g, dx[1] = %g, dx[2] = %g, dx[3] = %g\n", m4.vec_dx[0], m4.vec_dx[1], m4.vec_dx[2], m4.vec_dx[3]); 
	//} 
	//else 
	//{ 
	//  printf("Conjugate gradient method failed! f(x) = %g\n", fmin); 
	//} 
	//delete[]  xi; 
//}
}
