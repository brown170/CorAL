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
#include <iomanip>
#include <string>
#include "message.h"
#include "parametermap.h"
#include "cheezyparser.h"
#include "kernel_chooser.h"
#include "kernel.h"

using namespace std;

//-------------------------------------------------------------------------
// Main routines
//-------------------------------------------------------------------------

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"Simple code to pregenerate the kernel, expanded in Legendre polynomials, \n";
    cout<<"For use in Harmonic expansions of the sources and correlations.\n";
    cout<<"\nUsage: shark <option> [inputFile.dat]"<<endl;
    cout<<"    options -h, -help, --help all print this message, then quit"<<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** Simple HARmonic Kernels (SHARK) using CorAL ***"<<endl;
    cout << endl;
    
    bool got_file = false;

    MESSAGE << CMessage::warning;

    // Parse command line
    if (argc==1) getHelp();
    string paramFile("");
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="-help") getHelp();
        if (sarg=="--help") getHelp();
        if (sarg=="-h") getHelp();
        if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
        else {
            paramFile = sarg;
            got_file = true;
        }
    }    
    
    // Read in the input parameters from the argument on the command line
    if (!got_file) {
        MESSAGE<<"No inputFile parameter file given!!"<<ENDM_WARN;
        getHelp();
    }
    parameterMap inMap;
    parameter::ReadParsFromFile(inMap, paramFile);

    string kernel_particle1 = parameter::getS(inMap,"kernel_particle1","");
    string kernel_particle2 = parameter::getS(inMap,"kernel_particle2","");
    if ( inMap.hasKey("read_cache") ) MESSAGE << paramFile << " has read_cache key, I'm setting it to false"<<ENDM_WARN;
    parameter::set( inMap, "read_cache", false );
    if ( inMap.hasKey("use_cache") ) MESSAGE << paramFile << " has use_cache key, I'm setting it to true"<<ENDM_WARN;
    parameter::set( inMap, "use_cache", true );
    parameter::set( inMap, "param_filename", paramFile );
    

    cout << "Kernel generation settings:"<<inMap <<endl<<endl;
    
    chooseKernel( kernel_particle1, kernel_particle2, inMap );

    return true;
}
    
