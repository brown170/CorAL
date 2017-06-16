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
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "sf.h"
#include "constants.h"
using namespace std;


vector<double> bessj_zeros(int l, int nmax){
    vector<double> zeros;
    const double step_size(PI/1000.);
    double xlo(0.01),xhi(xlo+step_size);
    double ylo,yhi;
    while (zeros.size()<nmax) {
        ylo=jn(l,xlo);
        yhi=jn(l,xhi);
        if (ylo*yhi<0.0) zeros.push_back((xhi+xlo)/2.);
        xlo=xhi;
        xhi+=step_size;
    }
    return zeros;
}

int main(void){
    for (int l=0;l<4;l+=2){
        vector<double> zeros(bessj_zeros(l,10));
        cout << "l="<<l<<endl;
        for (unsigned int i=0;i<zeros.size();++i){
            cout << i << " " << zeros[i] << endl;
        }
    }
    return true;
}
