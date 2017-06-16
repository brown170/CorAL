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
#if !defined UTILS_H
#define      UTILS_H

#include <string>


//-------------------------------------------------------------------------
// Commonly used macros
//-------------------------------------------------------------------------
#define IS_ODD(n)        ((n) & 1)
#define IS_EVEN(n)       (!(IS_ODD(n)))
#define NEGONE_TO_THE(n) (IS_ODD(n) ? -1 : 1)
#define SIGN(x)          ((x) >= 0.0 ? 1 : -1)
#define MAXF(x,y)        maxof<double>(x,y)
#define MINF(x,y)        minof<double>(x,y)
#define MAXI(x,y)        maxof<int>(x,y)
#define MINI(x,y)        minof<int>(x,y)

//-------------------------------------------------------------------------
// Declarations for utility function templates.
//-------------------------------------------------------------------------

template <class T>
int iround(T x)           // Always rounds half-integers up
{
	int i = static_cast<int>(x);

	if (x >= static_cast<T>(0.)) {if (x - static_cast<T>(i) >= 0.5) i++;}
	else						 {if (x - static_cast<T>(i) < -0.5) i--;}

	return i;
}

template <class T>
inline T sign(T x)            // Returns 1. if iround(x) even, -1. if odd
{
	return (iround(x) % 2) ? static_cast<T>(-1.) : static_cast<T>(1.);
}

template <class T>
inline T minof(T x1, T x2)             // Returns lesser of x1, x2
{
	return (x1 < x2) ? x1 : x2;
}

template <class T>
inline T maxof(T x1, T x2)             // Returns greater of x1, x2
{
	return (x1 > x2) ? x1 : x2;
}

//-------------------------------------------------------------------------
// Small utility functions
//-------------------------------------------------------------------------
bool file_exists( const char* file_name );
bool file_exists( const std::string& file_name );
double mod(double x, double y);

#endif
