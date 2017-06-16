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
#include <fstream>
#include <iomanip>
#include <string>
#include "yasper.h"
#include "message.h"
#include "parametermap.h"
#include "cheezyparser.h"
#include "kernel_chooser.h"
#include "sou1d_legendre.h"
#include "sou1d_laguerre.h"
#include "sou1d_hermite.h"
#include "sou1d_chebyshev.h"
#include "sou1d_histo.h"
#include "sou1d_gauss.h"
#include "sou3d_ylm.h"
#include "sou1d_bsplines.h"
#include "kernel.h"
#include "basisfunc_imager1d.h"
#include "uncoupled_imager3d.h"
#include "constants.h"
#include "corr1d_histo.h"
#include "corr3d_ylm.h"
#include "convolution.h"

using namespace std;

//-------------------------------------------------------------------------
// Function prototypes
//-------------------------------------------------------------------------
void getHelp(void);
int main(int argc, char* argv[]);

//-------------------------------------------------------------------------
// Function implementations
//-------------------------------------------------------------------------

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"\n";
    cout<<"\n";
    cout<<"\nUsage: converts2c <option> [inputFile.dat]"<<endl;
    cout<<"    options -h, -help, --help all print this message, then quit"<<endl;
    cout<<"    -legendre       inputFile.dat uses Legendre polynomial basis" <<endl;
    cout<<"    -bspline        inputFile.dat uses Basis Spline basis"   <<endl;
    cout<<"    -laguerre       inputFile.dat uses Laguerre function basis"  <<endl;
    cout<<"    -chebyshev      inputFile.dat uses Chebyshev polynomial basis" <<endl;
    cout<<"    -histogram      inputFile.dat uses Histogram basis" <<endl;
    cout<<"    -hermite        inputFile.dat uses Hermite function basis"   <<endl;
    cout<<"    -gaussian       inputFile.dat is a 1d Gaussian"   <<endl;
    cout<<"    -3d             inputFile.dat is a 3d source"   <<endl;
    cout<<"    -convolve       Use convoluter, not imager for this (uncertainty not included)"   <<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** Widget to Convert Sources to Correlations (ConvertS2C) using CorAL ***"<<endl;
    cout << endl;
    
    bool got_file = false;
    bool three_d_source=false;
    bool convolve_with_kernel=false;

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
        if (sarg=="-3d") three_d_source=true;
        else if (sarg=="-convolve") convolve_with_kernel=true;
        else if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
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

    if (three_d_source) {
        if (convolve_with_kernel) MESSAGE<<"Not written yet!!  You have to do it yourself term by term!!"<<ENDM_FATAL;
        CCorrFtn3dSphr answer;
        parameterMap outMap;
        parameterMap sourceInMap;
        parameter::ReadParsFromFile(sourceInMap, parameter::getS(inMap,"source_file","input_source.dat"));
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it)
        {    
            if      ( *it == "-legendre" )    { 
                CSourceFtn3dSphr< CSourceFtn1dLegendrePoly > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledLegendrePolyImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            }
            else if ( *it == "-bspline" )     { 
                CSourceFtn3dSphr< CSourceFtn1dBSpline > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledBasisSplineImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            } 
            else if ( *it == "-laguerre" )    { 
                CSourceFtn3dSphr< CSourceFtn1dLaguerrePoly > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledLaguerrePolyImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            }
            else if ( *it == "-chebyshev" )   { 
                CSourceFtn3dSphr< CSourceFtn1dChebyshevPoly > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledChebyshevPolyImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            }
            else if ( *it == "-histogram" )   { 
                CSourceFtn3dSphr< CSourceFtn1dHisto > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledHistoImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            }
            else if ( *it == "-hermite" )     { 
                CSourceFtn3dSphr< CSourceFtn1dHermitePoly > souIn; 
                souIn.Read( sourceInMap );
                souIn.readTerms();
                UncoupledHermitePolyImager3d imager;
                imager.Read( inMap );
                imager.convertSourceToCorrelation( souIn, answer, inMap );
            }
            else if ( *it == "-gaussian" )    { MESSAGE<<"Not available in 3d: '"<<*it<<"'"<<ENDM_WARN; }
            else  MESSAGE<<"Unknown mode: '"<<*it<<"'"<<ENDM_WARN;
        }    
        answer.Write(outMap);
        parameter::WriteParsToFile(outMap, parameter::getS(inMap,"correlation_file","output_correlation.dat"));
        answer.writeTerms();
    }
    else {
        CCorrFtn1dHisto answer;
        parameterMap outMap;
        yasper::ptr< CSourceFtnBase > pExpansion; 
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it)
        {    
            if      ( *it == "-legendre" )  { pExpansion = new CSourceFtn1dLegendrePoly; }
            else if ( *it == "-bspline" )   { pExpansion = new CSourceFtn1dBSpline; }
            else if ( *it == "-laguerre" )  { pExpansion = new CSourceFtn1dLaguerrePoly; }
            else if ( *it == "-chebyshev" ) { pExpansion = new CSourceFtn1dChebyshevPoly; }
            else if ( *it == "-histogram" ) { pExpansion = new CSourceFtn1dHisto; }
            else if ( *it == "-hermite" )   { pExpansion = new CSourceFtn1dHermitePoly; }
            else if ( *it == "-gaussian" )  { pExpansion = new CGaussianSource; }
            else  MESSAGE<<"Unknown mode: '"<<*it<<"'"<<ENDM_WARN;
        }
        parameterMap sourceInMap;
        parameter::ReadParsFromFile(sourceInMap, parameter::getS(inMap,"source_file","input_source.dat"));
        pExpansion->Read( sourceInMap );
        if (convolve_with_kernel) {
            CBasisFunctionExpansion1d* pJunk = dynamic_cast<CBasisFunctionExpansion1d*>(pExpansion.GetRawPointer());
            if (!pJunk) throw MESSAGE << "Dead pointer pJunk"<<ENDM_FATAL;
            answer = convolute( *pJunk, inMap );
        } else {
            CBasisFuncImager1d imager;
            imager.Read( inMap );
            imager.convertSourceToCorrelation( *pExpansion, answer, inMap );
        }
        answer.Write(outMap);
        parameter::WriteParsToFile(outMap, parameter::getS(inMap,"correlation_file","output_correlation.dat"));
    }
    return true;
}
    
