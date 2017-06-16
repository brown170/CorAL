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
#include <iostream>
#include <cstdlib>
#include <cmath>

#include "sf.h"
#include "misc.h"
#include ""

//
// This code was superceeded with the gslbesstester.cc code since we have chosen to 
// use the GSL Bessel function 
//
// DAB 11/5/2009
//

using namespace std;

int main(){
  int n;
  double z,ans;
  
 TRYNEW:
  cout << "What is z?\n";
  cin >> z;
  ans=Bessel::j0(z);
  cout << "j0(" << z << ")= " << ans << "=" << sin(z)/z << endl;
  ans=Bessel::j1(z);
  cout << "j1(" << z << ")= " << ans << "=" << -cos(z)/z+sin(z)/(z*z) << endl;


  cout << "What is n? (Choose n>0)";
  cin >> n;
  ans=Bessel::jn(n,z);
  cout << "j(" << n << "," << z << ")= " << ans << endl;

  cout << "Enter 0 to quit, other to continue\n";
  cin >> n;
  if(n!=0) goto TRYNEW;

  return 0;
  
}

