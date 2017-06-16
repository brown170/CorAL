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
#include <iomanip>
#include "utils.h"
#include "cheezyparser.h"
#include "parametermap.h"
#include "sou1d_legendre.h"
#include "sou1d_laguerre.h"
#include "sou1d_hermite.h"
#include "sou1d_chebyshev.h"
#include "sou1d_histo.h"
#include "sou1d_gauss.h"
#include "sou1d_bsplines.h"
#include "yasper.h"
#include "bspline_imager1d.h"
#include "basisfunc_imager1d.h"
#include "sou3d_ylm.h"
#include "corr3d_ylm.h"
#include "corr3d_histo.h"
#include "uncoupled_imager3d.h"

// ------------- getHelp -------------
void getHelp(void){
    cout<<"\nUsage:  diver [mode] <input file>"<<endl;
    cout<<"    To get this message, use \"-help\" "<<endl;
    cout<<"        \"-bspline\"    : Use Basis Spline basis for source"<<endl;
    cout<<"        \"-chebyshev\"  : Use Chebyshev polynomial basis for source"<<endl;
    cout<<"        \"-hermite\"    : Use Hermite polynomial basis for source"<<endl;
    cout<<"        \"-histogram\"  : Use Histogram basis for source"<<endl;
    cout<<"        \"-laguerre\"   : Use Laguerre polynomial basis for source"<<endl;
    cout<<"        \"-legendre\"   : Use Legendre polynomial basis for source"<<endl;
    cout<<"        \"-norestore\"  : Don't restore the correlation"<<endl;
    cout<<"        \"-restore\"    : Restore the correlation"<<endl;  
    cout<<"        \"-1d\"         : Problem is 1d"<<endl;
    cout<<"        \"-3d\"         : Problem is 3d"<<endl;
    exit(0);
}

// ------------- main -------------
int main( int argc,  char* argv[] ){
    cout << "*** Demonstration InVERter: a CorAL imaging code ***"<<endl;
    cout << endl;

    // Parse command line
    if (argc==1) getHelp();
    vector<string> modeList;
    string inFile;
    bool got_file = false;
    bool restore  = true;
    bool three_d_problem = false;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if ( sarg == "-help" )  getHelp();
        if ( sarg == "-h" )     getHelp();
        if ( sarg == "--help" ) getHelp();
        if ( sarg == "-?" )     getHelp();
        if ( sarg == "-3d" ) three_d_problem = true;
        if ( sarg == "-norestore" ) restore=false;
        if ( sarg == "-restore" )   restore=true;
        if ( sarg == "-1d" ) three_d_problem = false;
        else if ( sarg.substr(0,1) == "-" ) modeList.push_back(sarg);
        else {
            inFile = sarg;
            got_file = true;
        }
    }
    if ( !got_file ) throw MESSAGE << "No imput file specified. Halting!!"<<ENDM_FATAL;
    if ( !file_exists( inFile ) ) throw MESSAGE << "File "<< inFile <<" not found. Halting!!"<<ENDM_FATAL;
   
    // Define & initialize the variables
    parameterMap settings, input_correlation_settings, imaged_source_settings, restored_correlation_settings;
    parameter::ReadParsFromFile(settings, inFile);
    string correlation_filename = parameter::getS( settings, "correlation_file" );
    string source_filename = parameter::getS( settings, "imaged_source_file" );
    string restored_correlation_filename = parameter::getS( settings, "restored_correlation_file" );
    parameter::ReadParsFromFile( input_correlation_settings, correlation_filename );

    if ( three_d_problem ) {
        // Define the input and output correlations
        cout << "Reading 3d correlation ..."<<endl;
        CCorrFtn3dSphr corrin(  "pi0", "pi0", false );
        CCorrFtn3dSphr corrout( "pi0", "pi0", false );
        cout << "    Input correlation: "<<correlation_filename<<endl;
        corrin.Read( input_correlation_settings );
        if (restore) corrout.Read( input_correlation_settings ); // don't worry, we'll overwrite the data
        corrin.readTerms();
        cout << endl;

        // Define the source and the appropriate imager
        yasper::ptr< CSourceFtnBase > pSource;
        yasper::ptr< CGeneralImager3d > pImager;
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it ) {
            if      ( *it == "-legendre" )  { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dLegendrePoly >; 
                pImager = new UncoupledLegendrePolyImager3d; 
            }
            else if ( *it == "-bspline" )   { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dBSpline >; 
                pImager = new UncoupledBasisSplineImager3d; 
            }
            else if ( *it == "-laguerre" )  { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dLaguerrePoly >; 
                pImager = new UncoupledLaguerrePolyImager3d; 
            }
            else if ( *it == "-chebyshev" ) { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dChebyshevPoly >; 
                pImager = new UncoupledChebyshevPolyImager3d; 
            }
            else if ( *it == "-histogram" ) { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dHisto >; 
                pImager = new UncoupledHistoImager3d; 
            }
            else if ( *it == "-hermite" )   { 
                pSource = new CSourceFtn3dSphr< CSourceFtn1dHermitePoly >; 
                pImager = new UncoupledHermitePolyImager3d; 
            }
            else if ( *it == "-3d" || *it == "-restore" || *it == "-norestore" ) {}
            else  MESSAGE<<"Unknown mode: '"<<*it<<ENDM_WARN;
        }
        // Set up the source and the imager
        cout << "Initializing 3d source ..."<<endl;
        pSource->Read( parameter::getMap( settings, "source_settings" ) );
    
        // Image & restore
        cout << "Imaging source ..."<<endl;
        pImager->convertCorrelationToSource( corrin, *pSource, settings );
        pSource->Write( imaged_source_settings );
        if (restore){
            cout << "Restoring correlation ..."<<endl;
            pImager->convertSourceToCorrelation( *pSource, corrout, settings );
            corrout.Write( restored_correlation_settings );
        }
     } 
    else {
        cout << "Reading 1d correlation ..."<<endl;
        // Define the input and output correlations
        CCorrFtn1dHisto corrin(  "pi0", "pi0", false );
        CCorrFtn1dHisto corrout( "pi0", "pi0", false );
        cout << "    Input correlation: "<<correlation_filename<<endl;
        corrin.Read( input_correlation_settings );
//        cout << input_correlation_settings <<endl;
        if (restore) corrout.Read( input_correlation_settings ); // don't worry, we'll overwrite the data

        // Define the source and the appropriate imager
        yasper::ptr< CSourceFtnBase > pSource;
        yasper::ptr< CBasisFuncImager1d > pImager;
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it ) {
            if      ( *it == "-legendre" )  { 
                pSource = new CSourceFtn1dLegendrePoly; 
                pImager = new CBasisFuncImager1d; 
            }
            else if ( *it == "-bspline" )   { 
                pSource = new CSourceFtn1dBSpline; 
                pImager = new CBasisSplineImager1d; 
            }
            else if ( *it == "-laguerre" )  { 
                pSource = new CSourceFtn1dLaguerrePoly; 
                pImager = new CBasisFuncImager1d; 
            }
            else if ( *it == "-chebyshev" ) { 
                pSource = new CSourceFtn1dChebyshevPoly; 
                pImager = new CBasisFuncImager1d; 
            }
            else if ( *it == "-histogram" ) { 
                pSource = new CSourceFtn1dHisto; 
                pImager = new CBasisFuncImager1d; 
            }
            else if ( *it == "-hermite" )   { 
                pSource = new CSourceFtn1dHermitePoly; 
                pImager = new CBasisFuncImager1d; 
            }
            else if ( *it == "-1d" || *it == "-restore" || *it == "-norestore" ) {}
            else  MESSAGE<<"Unknown mode: '"<<*it<<ENDM_WARN;
        }

        // Set up the source and the imager
        cout << "Initializing 1d source ..."<<endl;
        pSource->Read( parameter::getMap( settings, "source_settings" ) );
    
        // Image & restore
        cout << "Imaging source ..."<<endl;
        pImager->convertCorrelationToSource( corrin, *pSource, settings );
        pSource->Write( imaged_source_settings );

        if (restore) {
            cout << "Restoring correlation ..."<<endl;
            pImager->convertSourceToCorrelation( *pSource, corrout, settings );
            corrout.Write( restored_correlation_settings );
        }
    }
    
    // Write out the results to disk
    parameter::WriteParsToFile( imaged_source_settings, source_filename );
    if (restore) parameter::WriteParsToFile( restored_correlation_settings, restored_correlation_filename );
    
    return true;
}
