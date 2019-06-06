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
#include "sou1d_bsplines.h"
#include "sou1d_gauss.h"
#include "sou3d_ylm.h"
#include "corr3d_ylm.h"
#include "corr3d_histo.h"
#include "sou3d_histo.h"
#include "constants.h"

using namespace std;

//-------------------------------------------------------------------------
// Function prototypes
//-------------------------------------------------------------------------
void plotWidget1d( yasper::ptr< CObject1d > pExpansion, const parameterMap& inMap, double xmin, double xmax );
void plot3dSlices( yasper::ptr< CObject3d > pObj, const parameterMap& inMap );
void getHelp(void);
int main(int argc, char* argv[]);

//-------------------------------------------------------------------------
// Function implementations
//-------------------------------------------------------------------------

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"\n";
    cout<<"\nUsage: scplot <option> [inputFile.dat]"<<endl;
    cout<<"    options -h, -help, --help all print this message, then quit"<<endl;
    cout<<"    -legendre       Print out using Legendre polynomial basis" <<endl;
    cout<<"    -bspline        Print out using Basis Spline basis"   <<endl;
    cout<<"    -laguerre       Print out using Laguerre function basis"  <<endl;
    cout<<"    -chebyshev      Print out using Chebyshev polynomial basis" <<endl;
    cout<<"    -histogram      Print out using Histogram basis" <<endl;
    cout<<"    -hermite        Print out using Hermite function basis"   <<endl;
    cout<<"    -gaussian       Print out using using Gaussian class"   <<endl;
    cout<<"    -correlation    inputFile.dat is a correlation, only histogram basis available"   <<endl;
    cout<<"    -3dsphr         inputFile.dat is a 3d object, in spherical harmonics"<<endl;
    cout<<"    -3dcarthisto    inputFile.dat is a 3d Cartesian histogram object"<<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** Widget to PLOT Sources and Correlation (SCPLOT) using CorAL ***"<<endl;
    cout << "            (we need a catchy name for this code) "<<endl;
    cout << endl;

    bool got_file = false;
    bool is_3dsphr_object = false;
    bool is_3dcarthisto_object = false;

    // Parse command line
    if (argc==1) getHelp();
    string paramFile("");
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="-help") getHelp();
        if (sarg=="--help") getHelp();
        if (sarg=="-h") getHelp();
        if (sarg=="-3dsphr") is_3dsphr_object=true;
        else if (sarg=="-3dcarthisto") is_3dcarthisto_object=true;
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
    parameter::ReadParsFromFile( inMap, paramFile );

    if ( !is_3dsphr_object && !is_3dcarthisto_object ) {
        yasper::ptr< CObject1d > pExpansion;
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it)
        {
            if      ( *it == "-legendre" )    { pExpansion = new CLegendrePolynomialExpansion1d; }
            else if ( *it == "-bspline" )     { pExpansion = new CBasisSpline1d; }
            else if ( *it == "-laguerre" )    { pExpansion = new CLaguerrePolynomialExpansion1d; }
            else if ( *it == "-chebyshev" )   { pExpansion = new CChebyshevPolynomialExpansion1d; }
            else if ( *it == "-histogram" )   { pExpansion = new CHistogram1d; }
            else if ( *it == "-hermite" )     { pExpansion = new CHermiteFunctionExpansion1d; }
            else if ( *it == "-gaussian" )    { pExpansion = new CGaussianSource; }
            else if ( *it == "-correlation" ) { pExpansion = new CHistogram1d; }
            else  MESSAGE<<"Unknown mode: '"<<*it<<ENDM_WARN;
        }
        pExpansion->Read( inMap );
        plotWidget1d( pExpansion, inMap, 0.0, 100.0 );
    }
    else if ( is_3dsphr_object ){
        yasper::ptr< CObject3d > pObj;
        for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it)
        {
            if      ( *it == "-legendre" )    { pObj = new CSphericalHarmonicExpansion< CLegendrePolynomialExpansion1d >; }
            else if ( *it == "-bspline" )     { pObj = new CSphericalHarmonicExpansion< CBasisSpline1d >; }
            else if ( *it == "-laguerre" )    { pObj = new CSphericalHarmonicExpansion< CLaguerrePolynomialExpansion1d >; }
            else if ( *it == "-chebyshev" )   { pObj = new CSphericalHarmonicExpansion< CChebyshevPolynomialExpansion1d >; }
            else if ( *it == "-histogram" )   { pObj = new CSphericalHarmonicExpansion< CHistogram1d >; }
            else if ( *it == "-hermite" )     { pObj = new CSphericalHarmonicExpansion< CHermiteFunctionExpansion1d >; }
            else if ( *it == "-gaussian" )    { pObj = new CSphericalHarmonicExpansion< CGaussianSource >; }
            else if ( *it == "-correlation" ) { pObj = new CSphericalHarmonicExpansion< CHistogram1d >; }
            else  MESSAGE<<"Unknown mode: '"<<*it<<ENDM_WARN;
        }
        pObj->Read( inMap );
        pObj->readTerms();
        plot3dSlices( pObj, inMap );
    }
    else if ( is_3dcarthisto_object ){
        yasper::ptr< CObject3d > pObj = new CHistogram3d;
        pObj->Read( inMap );
        plot3dSlices( pObj, inMap );
    }
    return 0;
}

// ------------- plotWidget1d -------------
void plotWidget1d( yasper::ptr< CObject1d > pExpansion, const parameterMap& inMap, double xmin, double xmax ){
    string outFile = parameter::getS(inMap,"plot_file","output_plot.dat");
    cout << "    Plotting object to " << outFile << endl;
    double x;
    ofstream fout(outFile.c_str());
    for (int n=0;n<100;++n) {
        x = (xmax-xmin)*double(n)/100.0 + xmin;
        fout << x << " " << pExpansion->getValue(x)<<" "<<pExpansion->getError(x)<<endl;
    }
}

// ------------- plot3dSlices -------------
void plot3dSlices( yasper::ptr< CObject3d > pObj, const parameterMap& inMap ){
    string xFileName = parameter::getS(inMap,"plot_file_x","output_plot_x.dat");
    string yFileName = parameter::getS(inMap,"plot_file_y","output_plot_y.dat");
    string zFileName = parameter::getS(inMap,"plot_file_z","output_plot_z.dat");
    cout << "Writing X slice to plot file "<<xFileName<<endl;
    ofstream dataX( xFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataX << x << "  " << pObj->getValueCart(x,0.,0.) << "  " << pObj->getErrorCart(x,0.,0.)<<endl;
    }
    cout << "Writing Y slice to plot file "<<yFileName<<endl;
    ofstream dataY( yFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataY << x << "  " << pObj->getValueCart(0.,x,0.) << "  " << pObj->getErrorCart(0.,x,0.)<<endl;
    }
    cout << "Writing Z slice to plot file "<<zFileName<<endl;
    ofstream dataZ( zFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataZ << x << "  " << pObj->getValueCart(0.,0.,x) << "  " << pObj->getErrorCart(0.,0.,x)<<endl;
    }
}
