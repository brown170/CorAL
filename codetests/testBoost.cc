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
#include "misc.h"
#include <iostream>

int main( void ) {
    double u[4]={1.,0,0,.2};
    double gamma = 1.0/sqrt(1.0-u[1]*u[1]-u[2]*u[2]-u[3]*u[3]);
    for (int i=0;i<4;++i){u[i]*=gamma;}
    double p[4]={1.,0.,0.,0.};
    double pp[4];
    double gdiag[4]={1.,-1.,-1.,-1.};
    Misc::lorentz(u,p,pp);

    double u2 = 0.0;
    double p2 = 0.0;
    double pp2 = 0.0;

    for (int i=0;i<4;++i){
        u2+=u[i]*u[i]*gdiag[i];
        p2+=p[i]*p[i]*gdiag[i];
        pp2+=pp[i]*pp[i]*gdiag[i];
    }

    cout << "gamma: "<<gamma<<endl;
    cout << "u: ("; for (int i=0;i<4;++i) cout << u[i]<<",";cout << "),  u2: "<<u2<<endl;
    cout << "p: ("; for (int i=0;i<4;++i) cout << p[i]<<",";cout << "),  p2: "<<p2<<endl;
    cout << "p': ("; for (int i=0;i<4;++i) cout << pp[i]<<",";cout << "),  p'2: "<<pp2<<endl;

    return 1;
}
