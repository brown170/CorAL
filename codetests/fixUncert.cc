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
#import "cheezyparser.h"
#import "parametermap.h"
#import "linalg.h"
#import <iostream>

// compile with:
//    g++ fixUncert.cc -L../../coral_msu_repo/trunk/lib/ -lm -lgsl -lcoral -lcoralutils -I../../coral_msu_repo/trunk/include

int main(void){
    string infile, outfile;
    cout << "input file:"<<flush;
    cin >> infile;
    cout << endl;
    parameterMap m;
    parameter::ReadParsFromFile( m, infile );
    TNT::Array2D<double> covmtx;
    covmtx = stl2tntMat(parameter::getM(m,"covmtx_datablock"));
    parameter::set(m,"covmtx_datablock",tnt2stlMat(10000.0*covmtx));
    //cout << m;
    cout << "output file:"<<flush;
    cin >> outfile;
    cout << endl;
    parameter::WriteParsToFile( m, outfile );
    return true;
}
