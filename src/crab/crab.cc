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
#include "parametermap.h"
#include "cheezyparser.h"
#include "oscar.h"
#include "oscar_source_generator_1d.h"
#include "oscar_correlation_generator_1d.h"
#include "oscar_source_generator_3dcart.h"
#include "oscar_source_generator_3dsphr.h"
#include "oscar_correlation_generator_3dcart.h"
#include "message.h"
#include "basisfunc_imager1d.h"
#include "expander.h"
#include "uncoupled_imager3d.h"
#include "corr3d_ylm.h"
#include "tnt_array1d.h"
#include <iostream>
#include <iomanip>

//-------------------------------------------------------------------------
// Function declarations
//-------------------------------------------------------------------------

// ------------- wrappers around source generators -------------
CSourceFtn1dHisto getSource1d( vector<COSCARLine> particleList, parameterMap m );
CSourceFtn1dHisto getSource1dOldWay( vector<COSCARLine> particleList, parameterMap m );
CSourceFtn3dHisto getSource3dCart( vector<COSCARLine> particleList, parameterMap m );
CSourceFtn3dSphr< CSourceFtn1dHisto >
                  getSource3dSphr( vector<COSCARLine> particleList, parameterMap m );

// ------------- wrappers around correlation generators -------------
CCorrFtn1dHisto   getCorrelation1d( vector<COSCARLine> particleList, parameterMap m );
CCorrFtn1dHisto   getCorrelation1dOldWay( const CSourceFtn1dHisto& souin, parameterMap m );
CCorrFtn3dHisto   getCorrelation3dCart( vector<COSCARLine> particleList, parameterMap m );
CCorrFtn3dSphr    getCorrelation3dSphr( vector<COSCARLine> particleList, parameterMap m );

// ------------- slice plotting routines -------------
void plot3dSlices( CObject3d& obj, string xFileName, string yFileName, string zFileName );
void plot3dSlices( CHistogram3d& obj, string xFileName, string yFileName, string zFileName );
void plot1d( CObject1d& obj, string outFileName );
void plot1d( CHistogram1d& obj, string outFileName );

// ------------- driver routines -------------
void make1dSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots );
void make1dSourceAndCorrelationOldWay( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots );
void make3dCartSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots );
void make3dCartExpandedSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots );
void make3dSphrSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots );
void makeCRABCorrelations( vector<COSCARLine> particleList, parameterMap& inMap, bool make_plots );
void makeSourceMakerSources( vector<COSCARLine> particleList, parameterMap& inMap, bool make_plots );

//-------------------------------------------------------------------------
// Main routines
//-------------------------------------------------------------------------

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"\nUsage: crab <mode> [inputFile.dat]"<<endl;
    cout<<"    Valid modes are: -1d           : Generate 1d source and correlation"<<endl;
    cout<<"                     -1doldway     : Generate 1d source and correlation, using old code for source"<<endl;
    cout<<"                     -3dcart       : Generate 3d source and correlation directly in Cartesian"<<endl;
    cout<<"                     -3dcart_ex    : Generate 3d source and correlation directly in Cartesian, then expand"<<endl;
    cout<<"                     -3dsphr       : Generate 3d source and correlation directly in Ylm's"<<endl;
    cout<<"                     -crab         : Generate 1d & 3d correlations and output them in the style of Scott Pratt's CRAB code"<<endl;
    cout<<"                     -sourcemaker  : Generate 1d source & 3d source slices and output them in the style of Dave Brown's SourceMaker code"<<endl;
    cout<<"                     -skip_correlation  : Skip generation of the correlations if accompanying modes would generate them"<<endl;
    cout<<"                     -make_plots   : Make plots for all the sources and correlations generated"<<endl;
    cout<<"    -h, -help, --help              : Prints this message, then quits"<<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** sample 1d implementation of the CoRrelation AfterBurner (CRAB) using CorAL ***"<<endl;
    cout << endl;

    MESSAGE << CMessage::warning;
    //(MESSAGE).show_info=true;

    bool skip_correlation=false;
    bool make_plots=false;

    // Parse command line
    if (argc==1) getHelp();
    string paramFile("");
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="-help") getHelp();
        if (sarg=="--help") getHelp();
        if (sarg=="-h") getHelp();
        if (sarg=="-skip_correlation") {
            cout << "Will skip correlation function generation\n"<<endl;
            skip_correlation = true;
        }
        else if (sarg=="-make_plots") {
            cout << "Will make files suitable for plotting\n"<<endl;
            make_plots = true;
        }
        else if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
        else paramFile = sarg;
    }

    // Read in the input parameters from the argument on the command line
    if (paramFile=="") {
        MESSAGE<<"No inputFile parameter file given!!"<<ENDM_FATAL;
        getHelp();
    }
    parameterMap inMap;
    ReadParsFromFile(inMap, paramFile);

    // Read in the OSCAR formatted data
    string oscarFile = inMap.getItem( "oscar_file", (string)"thermal_gauss01.dat" );
    vector<COSCARLine> particleList = filterPID(
        readOSCARFile(oscarFile),
        inMap.getItem( "pid", 211 )
    );


    // Generate the correlation from the source, then put hte results in the map
    for (vector<string>::iterator mode=modeList.begin(); mode!=modeList.end(); ++mode) {
        if      (*mode == "-1d")                        make1dSourceAndCorrelation( particleList, inMap, skip_correlation, make_plots );
        else if (*mode == "-1doldway")                  make1dSourceAndCorrelationOldWay( particleList, inMap, skip_correlation, make_plots );
        else if (*mode == "-3dcart")                    make3dCartSourceAndCorrelation( particleList, inMap, skip_correlation, make_plots );
        else if (*mode == "-3dcart_ex")                 make3dCartExpandedSourceAndCorrelation( particleList, inMap, skip_correlation, make_plots  );
        else if (*mode == "-3dsphr")                    make3dSphrSourceAndCorrelation( particleList, inMap, skip_correlation, make_plots );
        else if (*mode == "-crab" && !skip_correlation) makeCRABCorrelations( particleList, inMap, make_plots );
        else if (*mode == "-sourcemaker")               makeSourceMakerSources( particleList, inMap, make_plots );
        else MESSAGE<< "Unknown mode: "<<*mode<<ENDM_WARN;
        cout << endl;
    }
    return 0;
}

//-------------------------------------------------------------------------
// Function definitions
//-------------------------------------------------------------------------

// ------------- make1dSourceAndCorrelation -------------
void make1dSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots ){
    cout << "\nGenerating 1d source"<<endl;
    CSourceFtn1dHisto souin = getSource1d( particleList, inMap );
    cout << endl;
    parameterMap coutMap, soutMap;
    souin.Write(soutMap);
    WriteParsToFile(soutMap,parameter::getS(inMap, "source_file", "output_source_1d.dat"));
    if (make_plots){
        plot1d( souin, getS( inMap, "source_plot_filename", "output_source_1d_plot.dat" ));
    }
    if (!skip_correlation) {
        cout << "\nGenerating 1d correlation from source"<<endl;
        CCorrFtn1dHisto answer = getCorrelation1d( particleList, inMap );
        cout << endl;
        answer.Write(coutMap);
        WriteParsToFile(coutMap,parameter::getS(inMap, "correlation_file", "output_correlation_1d.dat"));
        if (make_plots){
            plot1d( answer, getS( inMap, "correlation_plot_filename", "output_correlation_1d_plot.dat" ));
        }
    }
}

// ------------- make1dSourceAndCorrelationOldWay -------------
void make1dSourceAndCorrelationOldWay( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots ){
    cout << "\nGenerating 1d source, the old way"<<endl;
    parameterMap souMapIn;
    souMapIn = parameter::getMap(inMap, "source_settings");
    CSourceFtn1dHisto soualt = getSource1dOldWay( particleList, souMapIn );
    cout << endl;
    parameterMap coutMap, soutMap;
    soualt.Write(soutMap);
    WriteParsToFile(soutMap,parameter::getS(inMap, "source_file", "output_source_1d.dat"));
    if (make_plots){
        plot1d( soualt, getS( inMap, "source_plot_filename", "output_source_1d_plot.dat" ));
    }
    if (!skip_correlation) {
        cout << "\nGenerating 1d correlation from source"<<endl;
        CCorrFtn1dHisto answer;
        CBasisFuncImager1d imager;
        imager.convertSourceToCorrelation(soualt, answer, inMap);
        cout << endl;
        answer.Write(coutMap);
        WriteParsToFile(coutMap,parameter::getS(inMap, "correlation_file", "output_correlation_1d.dat"));
        if (make_plots){
            plot1d( answer, getS( inMap, "correlation_plot_filename", "output_correlation_1d_plot.dat" ));
        }
    }
}

// ------------- make3dCartSourceAndCorrelation -------------
void make3dCartSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots ){
    cout << "\nGenerating 3d source, in Cartesian"<<endl;
    CSourceFtn3dHisto sanswer = getSource3dCart( particleList, inMap );
    cout << endl;
    parameterMap coutMap, soutMap;
    sanswer.Write(soutMap);
    WriteParsToFile( soutMap, parameter::getS(inMap, "source_file", "output_source_3d.dat") );
    if (make_plots){
        plot3dSlices( sanswer,
            getS( inMap, "source_3d_out_plot_filename", "output_source_3d_out_plot.dat" ),
            getS( inMap, "source_3d_side_plot_filename", "output_source_3d_side_plot.dat" ),
            getS( inMap, "source_3d_long_plot_filename", "output_source_3d_long_plot.dat" )
        );
    }
    if (!skip_correlation) {
        cout << "\nGenerating 3d correlation, in Cartesian"<<endl;
        CCorrFtn3dHisto answer = getCorrelation3dCart( particleList, inMap );
        cout << endl;
        answer.Write(coutMap);
        WriteParsToFile( coutMap, parameter::getS(inMap, "correlation_file", "output_correlation_3d.dat") );
        if (make_plots){
            plot3dSlices( answer,
                getS( inMap, "correlation_3d_out_plot_filename", "output_correlation_3d_out_plot.dat" ),
                getS( inMap, "correlation_3d_side_plot_filename", "output_correlation_3d_side_plot.dat" ),
                getS( inMap, "correlation_3d_long_plot_filename", "output_correlation_3d_long_plot.dat" )
            );
        }
    }
}

// ------------- make3dCartExpandedSourceAndCorrelation -------------
void make3dCartExpandedSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots ){
    cout << "\nGenerating 3d source in Cartesian coordinates"<<endl;
    CSourceFtn3dHisto souin3d = getSource3dCart( particleList, inMap );
    cout << "\nExpanding 3d source into Ylm's"<<endl;
    CSourceFtn3dSphr< CSourceFtn1dHisto > soutemp = expand( souin3d, inMap );
    cout << endl;
    parameterMap coutMap, soutMap, souinMap;
    souin3d.Write(souinMap);
    soutemp.Write(soutMap);
    WriteParsToFile( souinMap, parameter::getS(inMap, "source_file_cart", "output_source_3d_cart.dat") );
    WriteParsToFile( soutMap, parameter::getS(inMap, "source_file", "output_source_3d_sphr.dat") );
    soutemp.writeTerms();
    if (make_plots){
        cout << "Do nothing for plotting"<<endl;
    }
    if (!skip_correlation) {
        cout << "\nGenerating 3d correlation from source"<<endl;
        CCorrFtn3dSphr answer;
        UncoupledHistoImager3d imager;
        imager.convertSourceToCorrelation(soutemp, answer, inMap);
        cout << endl;
        answer.Write(coutMap);
        WriteParsToFile( coutMap, parameter::getS(inMap, "correlation_file", "output_correlation_3d_cart.dat") );
        answer.writeTerms( );
        if (make_plots){
            cout << "Do nothing for plotting"<<endl;
        }
    }
}

// ------------- make3dSphrSourceAndCorrelation -------------
void make3dSphrSourceAndCorrelation( vector<COSCARLine> particleList, parameterMap& inMap, bool skip_correlation, bool make_plots ){
    cout << "\nGenerating 3d source directly in Ylm's"<<endl;
    CSourceFtn3dSphr< CSourceFtn1dHisto > soualt3d = getSource3dSphr( particleList, inMap );
    cout << endl;
    parameterMap soutMap, coutMap;
    soualt3d.Write(soutMap);
    soualt3d.writeTerms();
    WriteParsToFile( soutMap, getS( inMap, "source_file", "output_source_3dsphr.dat" ) );
    if (make_plots){
        plot3dSlices( soualt3d,
            getS( inMap, "source_3d_out_plot_filename", "output_source_3d_out_plot.dat" ),
            getS( inMap, "source_3d_side_plot_filename", "output_source_3d_side_plot.dat" ),
            getS( inMap, "source_3d_long_plot_filename", "output_source_3d_long_plot.dat" )
        );
    }
    if (!skip_correlation) {
        cout << "\nGenerating 3d correlation from source"<<endl;
        CCorrFtn3dSphr answer;
        UncoupledHistoImager3d imager;
        imager.convertSourceToCorrelation(soualt3d, answer, inMap);
        cout << endl;
        answer.Write(coutMap);
        answer.writeTerms();
        WriteParsToFile( coutMap, getS( inMap, "correlation_file", "output_correlation_3dsphr.dat" ) );
        if (make_plots){
            plot3dSlices( answer,
                getS( inMap, "correlation_3d_out_plot_filename", "output_correlation_3d_out_plot.dat" ),
                getS( inMap, "correlation_3d_side_plot_filename", "output_correlation_3d_side_plot.dat" ),
                getS( inMap, "correlation_3d_long_plot_filename", "output_correlation_3d_long_plot.dat" )
            );
        }
    }
}

// ------------- makeCRABCorrelations -------------
void makeCRABCorrelations( vector<COSCARLine> particleList, parameterMap& inMap, bool make_plots ){
    cout << "\nGenerating 1d correlation like CRAB used to do"<<endl;
    COSCARCorrelationGenerator1d corr1dgen;
    CCorrFtn1dHisto corr1d(corr1dgen.generateCorrelation(particleList, inMap));
    cout << endl;
    cout << "\nGenerating 3d correlation like CRAB used to do"<<endl;
    COSCARCorrelationGenerator3dCart corr3dgen;
    CCorrFtn3dHisto corr3d(corr3dgen.generateCorrelation(particleList, inMap));
    cout << endl;
    TNT::Array1D< int > numPairs1d = corr1dgen.pairCount;
    TNT::Array1D< int > numPairs3d = corr3dgen.pairCount;
    cout << "\nOutputting 1d and 3d correlations like CRAB used to do"<<endl;
    ofstream outFile1d(parameter::getS(inMap, "crab_correlation_1d_filename", "correlation_qinv.dat").c_str());
    outFile1d << "!As a function of reduced momentum, k=(p2-p1)/2:\n"
         << "!k(MeV/c)   C(k)     +/- numcounts&dencounts"<<endl;
    for (int i=0; i<corr1d.ndata; ++i)
    {
        outFile1d<< setw(6) << setprecision(3) << corr1d.midBin(i) << "    "
            << setw(6) << setprecision(3) << corr1d.data[i]   << "  "
            << setw(9) << setprecision(6) << corr1d.uncert[i] << "  "
            << numPairs1d[i]    << "  " << 0 << endl;
    }
    ofstream outFile3d(parameter::getS(inMap, "crab_correlation_3d_filename", "correlation_qinv3d.dat").c_str());
    outFile3d << "!As a function of reduced momentum, k=(p2-p1)/2:\n"
         << "! kx     ky     kz        C(k)    +/-  numcounts&dencounts"<<endl;
    for (int ix=0; ix<corr3d.nx; ++ix)
        for (int iy=0; iy<corr3d.ny; ++iy)
            for (int iz=0; iz<corr3d.nz; ++iz)
            {
                int i = corr3d.whatIndex(ix,iy,iz);
                outFile3d
                    << setw(6) << setprecision(3) << corr3d.midBinX(ix) << "    "
                    << setw(6) << setprecision(3) << corr3d.midBinY(iy) << "    "
                    << setw(6) << setprecision(3) << corr3d.midBinZ(iz) << "    "
                    << setw(6) << setprecision(3) << corr3d.data[i]     << "  "
                    << setw(9) << setprecision(6) << corr3d.uncert[i]   << "  "
                    << numPairs3d[i]    << "  " << 0 << endl;
            }
    if (make_plots) {
        plot1d( corr1d, getS( inMap, "correlation_plot_filename", "output_correlation_1d_plot.dat" ) );
        plot3dSlices( corr3d,
            getS( inMap, "correlation_3d_out_plot_filename", "output_correlation_3d_out_plot.dat" ),
            getS( inMap, "correlation_3d_side_plot_filename", "output_correlation_3d_side_plot.dat" ),
            getS( inMap, "correlation_3d_long_plot_filename", "output_correlation_3d_long_plot.dat" )
        );
    }
}

// ------------- makeSourceMakerSources -------------
void makeSourceMakerSources( vector<COSCARLine> particleList, parameterMap& inMap, bool make_plots ){
    cout << "\nGenerating 1d source slices like SourceMaker used to do"<<endl;
    CSourceFtn3dHisto source3d = getSource3dCart( particleList, inMap );
    cout << endl;
    cout << "\nGenerating 3d source slices like SourceMaker used to do"<<endl;
    CSourceFtn1dHisto source1d = getSource1d( particleList, inMap );
    if (make_plots) {
        plot1d( source1d, getS( inMap, "source_1d_plot_filename", "output_source_1d_plot.dat" ));
        cout << endl;
        plot3dSlices( source3d,
            getS( inMap, "source_3d_out_plot_filename", "output_source_3d_out_plot.dat" ),
            getS( inMap, "source_3d_side_plot_filename", "output_source_3d_side_plot.dat" ),
            getS( inMap, "source_3d_long_plot_filename", "output_source_3d_long_plot.dat" ) );
    }
}

// ------------- getSource1dOldWay -------------
//! Original source generator
CSourceFtn1dHisto getSource1dOldWay( vector<COSCARLine> particleList, parameterMap m ){
    CSourceFtn1dHisto souin(getOSCARSource1d(particleList, m));
    return souin;
}

// ------------- getSource1d -------------
//! New source generator for testing
CSourceFtn1dHisto getSource1d( vector<COSCARLine> particleList, parameterMap m ){
    COSCARSourceGenerator1d sougen;
    CSourceFtn1dHisto soualt(sougen.generateSource(particleList, m));
    return soualt;
}

// ------------- getCorrelation1d -------------
//! New correlation generator for testing
CCorrFtn1dHisto getCorrelation1d( vector<COSCARLine> particleList, parameterMap m ){
    COSCARCorrelationGenerator1d corrgen;
    CCorrFtn1dHisto corralt(corrgen.generateCorrelation(particleList, m));
    return corralt;
}

// ------------- getCorrelation1dOldWay -------------
//! Old correlation generator for testing, uses imager
CCorrFtn1dHisto getCorrelation1dOldWay(  const CSourceFtn1dHisto& souin, parameterMap m ){
    CBasisFuncImager1d imager;
    CCorrFtn1dHisto answer;
    imager.convertSourceToCorrelation( souin, answer, m );
    return answer;
}

// ------------- getSource3dSphr -------------
CSourceFtn3dSphr< CSourceFtn1dHisto > getSource3dSphr( vector<COSCARLine> particleList, parameterMap m ){
    COSCARSourceGenerator3dSphr sougen;
    CSourceFtn3dSphr< CSourceFtn1dHisto > souin(sougen.generateSource(particleList, m));
    return souin;
}

// ------------- getSource3dCart -------------
CSourceFtn3dHisto getSource3dCart( vector<COSCARLine> particleList, parameterMap m ){
    COSCARSourceGenerator3dCart sougen;
    CSourceFtn3dHisto souin(sougen.generateSource(particleList, m));
    return souin;
}

// ------------- getCorrelation3dCart -------------
CCorrFtn3dHisto getCorrelation3dCart( vector<COSCARLine> particleList, parameterMap m ){
    COSCARCorrelationGenerator3dCart sougen;
    CCorrFtn3dHisto souin(sougen.generateCorrelation(particleList, m));
    return souin;
}

// ------------- plot1d -------------
void plot1d( CObject1d& obj, string outFileName )
{
    cout << "Writing 1d plot file "<<outFileName<<endl;
    ofstream dataOut( outFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataOut << x << "  " << obj.getValue(x) << "  " << obj.getError(x)<<endl;
    }
}

// ------------- plot1d -------------
void plot1d( CHistogram1d& obj, string outFileName )
{
    cout << "Writing 1d plot file "<<outFileName<<endl;
    ofstream dataOut( outFileName.c_str() );
    for (int i=0; i<obj.ndata; ++i){
        double x = obj.midBin(i);
        dataOut << x << "  " << obj.getValue(x) << "  " << obj.getError(x)<<endl;
    }
}

// ------------- plot3dSlices -------------
void plot3dSlices( CObject3d& obj, string xFileName, string yFileName, string zFileName ){
    cout << "Writing X slice to plot file "<<xFileName<<endl;
    ofstream dataX( xFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataX << x << "  " << obj.getValueCart(x,0.,0.) << "  " << obj.getErrorCart(x,0.,0.)<<endl;
    }
    cout << "Writing Y slice to plot file "<<yFileName<<endl;
    ofstream dataY( yFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataY << x << "  " << obj.getValueCart(0.,x,0.) << "  " << obj.getErrorCart(0.,x,0.)<<endl;
    }
    cout << "Writing Z slice to plot file "<<zFileName<<endl;
    ofstream dataZ( zFileName.c_str() );
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataZ << x << "  " << obj.getValueCart(0.,0.,x) << "  " << obj.getErrorCart(0.,0.,x)<<endl;
    }
}

// ------------- plot3dSlices -------------
void plot3dSlices( CHistogram3d& obj, string xFileName, string yFileName, string zFileName ){
    cout << "Writing X slice to plot file "<<xFileName<<endl;
    ofstream dataX( xFileName.c_str() );
    for (int i=0; i<obj.nx; ++i){
        double x = obj.midBinX(i);
        dataX << x << "  " << obj.getValueCart(x,0.,0.) << "  " << obj.getErrorCart(x,0.,0.)<<endl;
    }
    cout << "Writing Y slice to plot file "<<yFileName<<endl;
    ofstream dataY( yFileName.c_str() );
    for (int i=0; i<obj.ny; ++i){
        double x = obj.midBinY(i);
        dataY << x << "  " << obj.getValueCart(0.,x,0.) << "  " << obj.getErrorCart(0.,x,0.)<<endl;
    }
    cout << "Writing Z slice to plot file "<<zFileName<<endl;
    ofstream dataZ( zFileName.c_str() );
    for (int i=0; i<obj.nz; ++i){
        double x = obj.midBinZ(i);
        dataZ << x << "  " << obj.getValueCart(0.,0.,x) << "  " << obj.getErrorCart(0.,0.,x)<<endl;
    }
}
