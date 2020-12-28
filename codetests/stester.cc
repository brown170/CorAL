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
#include "utils.h"
#include "cheezyparser.h"
#include "parametermap.h"
#include "sou1d_gauss.h"
#include "sou3d_ylm.h"
#include "corr3d_ylm.h"
#include "corr3d_histo.h"
#include "histogram1d.h"
#include "expander.h"
#include "convolution.h"
#include "imager.h"
#include "sou1d_bsplines.h"
#include "sou1d_histo.h"
#include "oscar.h"
#include <fstream>
#include <iostream>
#include <iomanip>

bool parameterMap_tester(void);
bool object_io_tester(void);
bool harmonic_tester(void);
bool imager_tester(void);
bool legendre_imager_tester(void);
bool convolution_tester(void);
bool legendre_3dimager_tester(void);

// ------------- getHelp -------------
void getHelp(void){
    cout<<"\nUsage: stester <test>"<<endl;
    cout<<"    To get this message, use \"-help\" "<<endl;
    cout<<"    Valid tests are: -pmap       : Test read/write/assignment of parameterMap"<<endl;
    cout<<"                     -objio      : Read/write various objects"<<endl;
    cout<<"                     -harm       : Initialize a 3d correlation represented in a Spherical Harmonic basis"<<endl;
    cout<<"                     -image      : Simple basis spline imaging test in 1d"<<endl;
    cout<<"                     -legimage   : Simple Legendre polynomial imaging test in 1d"<<endl;
    cout<<"                     -legimage3d : Legendre polynomial imaging test in 3d"<<endl;
    cout<<"                     -conv       : Convolute a Gaussian with a kernel"<<endl;
    exit(0);
}

// ------------- main -------------
int main( int argc,  char* argv[] ){
    cout << "*** simple CorAL regression tests ***"<<endl;
    cout << endl;

    // Parse command line
    if (argc==1) getHelp();
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="-help") getHelp();
        if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
    }    
    
    // Generate the correlation from the source, then put hte results in the map
    for (vector<string>::iterator mode=modeList.begin(); mode!=modeList.end(); ++mode) {
        if ( *mode == "-pmap" )       return parameterMap_tester();
        if ( *mode == "-objio" )      return object_io_tester();
        if ( *mode == "-harm" )       return harmonic_tester();
        if ( *mode == "-image" )      return imager_tester();
        if ( *mode == "-legimage" )   return legendre_imager_tester();
        if ( *mode == "-legimage3d" ) return legendre_3dimager_tester();
        if ( *mode == "-conv" )       return convolution_tester();
    }
    return false;
}



// ------------- convolution_tester -------------
bool convolution_tester(void){
    parameterMap m_in;
    parameterMap m_sou_in;
    parameterMap m_conv_in;
    parameterMap m_corr_in;
    parameterMap m_kern_in;
    // set up input source
    parameter::set(m_sou_in,"particle1",string("pi+"));
    parameter::set(m_sou_in,"particle2",string("pi+"));
    parameter::set(m_sou_in,"width",5.0);
    parameter::set(m_sou_in,"height",1.0);
    parameter::set(m_sou_in,"l",0);
    parameter::set(m_sou_in,"m",0);
    parameter::set(m_sou_in,"realpart",true);
    // set up for the convoluter
    parameter::set(m_conv_in,"rmin",0.0);
    parameter::set(m_conv_in,"rmax",50.0);
    // set up for the kernel
    parameter::set(m_kern_in,"hbt_only",false);
    parameter::set(m_kern_in,"kernel_cache",string("kdata/pi+pi+/"));
    parameter::set(m_kern_in,"use_cache",true);
    parameter::set(m_kern_in,"read_cache",true);
    parameter::set(m_kern_in,"param_filename",string("test_data/conv/wfparameters_new.dat"));
    // set up for the correlation we're creating
    parameter::set(m_corr_in,"fixed_width_bins",true);
    parameter::set(m_corr_in,"dx",3.0);
    parameter::set(m_corr_in,"xoffset",1.50);
    parameter::set(m_corr_in,"ndata",40);
    // load 'em into the main map
    parameter::set(m_in,"convoluter_settings",m_conv_in);
    parameter::set(m_in,"correlation_settings",m_corr_in);
    parameter::set(m_in,"kernel_settings",m_kern_in);

    CGaussianSource source;
    source.Read(m_sou_in);
    parameterMap m_test;
    source.Write(m_test);
    cout << "_____Gaussian Source Parameters_____";
    cout << m_test <<endl;
    cout << "source.getValue(2.0)="<<source.getValue(2.0)<<endl <<endl;
    
    CCorrFtn1dHisto correlation(convolute(source,m_in));
    for (int i=0;i<correlation.ndata;++i){
        correlation.uncert[i]=correlation.data[i]*0.1;
    }
    
    parameterMap m_out;
    correlation.Write(m_out);
    ofstream corrout("pion_correlation.dat");
    corrout << m_out;

    cout << "\nSample points from correlation:\nq   C(q)"<<endl;
    for (int i=0;i<10;++i){
        double q=5.0*(double(i)+0.5);
        cout << q << " " << correlation.getValue(q) << endl;
    }
    
    return true;
}



// ------------- imager_tester -------------
bool imager_tester(void){
    if ( !file_exists( "test_data/image/pion_correlation.dat") ) 
        throw MESSAGE << "File "<< "pion_correlation.dat" <<" not found. Halting!!"<<ENDM_FATAL;
    CCorrFtn1dHisto corrin("pi0","pi0",false);
    parameterMap corrin_params;
    ifstream corrin_file( "test_data/image/pion_correlation.dat" );
    corrin_file >> corrin_params;
    cout << "_________Input Correlation_________\n"<<corrin_params<<endl;
    corrin.Read(corrin_params);

    if ( !file_exists( "test_data/image/image.com") ) 
        throw MESSAGE << "File "<< "image.com" <<" not found. Halting!!"<<ENDM_FATAL;
    parameterMap image_com;
    ifstream image_file("test_data/image/image.com");
    image_file>>image_com;
    cout << "_________Imaging Controls_________\n" << image_com;
    
    cout << "_________Imaging_________\n";
    CSourceFtn1dBSpline souout(image(corrin,image_com));
    
    parameterMap sou_map;
    souout.Write(sou_map);
    cout << "_________Imaged Source_________\n" << sou_map;
    
    ofstream sourceOut("pion_source.dat");
    sourceOut << sou_map;
    
    ofstream dataOut("data.dat");
    for (int i=0; i<32; ++i){
        double x = (double)i;
        dataOut << x << "  " << souout.getValue(x) << "  " << souout.getError(x)<<endl;
    }
    
    return true;
}

// ------------- legendre_imager_tester -------------
bool legendre_imager_tester(void){
    if ( !file_exists( "test_data/legimage/pion_correlation.dat") ) 
        throw MESSAGE << "File "<< "test_data/legimage/pion_correlation.dat" <<" not found. Halting!!"<<ENDM_FATAL;
    CCorrFtn1dHisto corrin("pi0","pi0",false);
    parameterMap corrin_params;
    ifstream corrin_file( "test_data/legimage/pion_correlation.dat" );
    corrin_file >> corrin_params;
    // cout << "_________Input Correlation_________\n"<<corrin_params<<endl;
    corrin.Read(corrin_params);

    if ( !file_exists( "test_data/legimage/image.com") ) 
        throw MESSAGE << "File "<< "test_data/legimage/image.com" <<" not found. Halting!!"<<ENDM_FATAL;
    parameterMap image_com;
    ifstream image_file("test_data/legimage/image.com");
    image_file>>image_com;
    //cout << "_________Imaging Controls_________\n" << image_com;
    
    cout << "_________Imaging_________\n";
    CSourceFtn1dLegendrePoly souout(image_le(corrin,image_com));
    
    parameterMap sou_map;
    souout.Write(sou_map);
    //cout << "_________Imaged Source_________\n" << sou_map;
    
    ofstream sourceOut("pion_legendre_source.dat");
    sourceOut << sou_map;
    
    ofstream dataOut("legendre_data.dat");
    for (int i=0; i<52; ++i){
        double x = (double)i;
        dataOut << x << "  " << souout.getValue(x) << "  " << souout.getError(x)<<endl;
    }
    
    return true;
}

// ------------- legendre_3dimager_tester -------------
bool legendre_3dimager_tester(void){
    parameterMap com, image_com;
//    ifstream f("HBTOnlyTestData_corr.coral");
    ifstream f("baselineModelData_corr.coral");
    f >> com;
    cout << com;

    CCorrFtn3dSphr c;
    c.Read(com);
    c.readTerms();

    ifstream image_file("image.com");
    image_file>>image_com;
    cout << image_com;
    
    CSourceFtn3dSphr<CSourceFtn1dLegendrePoly> s = image_l(c,image_com);
    
    ofstream dataOutO("3d_legendre_data_O.dat");
    for (int i=0; i<62; ++i){
        double x = (double)i;
        dataOutO << x << "  " << s.getValueCart(x,0.,0.) << "  " << s.getErrorCart(x,0.,0.)<<endl;
    }
    ofstream dataOutS("3d_legendre_data_S.dat");
    for (int i=0; i<62; ++i){
        double x = (double)i;
        dataOutS << x << "  " << s.getValueCart(0.,x,0.) << "  " << s.getErrorCart(0.,x,0.)<<endl;
    }
    ofstream dataOutL("3d_legendre_data_L.dat");
    for (int i=0; i<62; ++i){
        double x = (double)i;
        dataOutL << x << "  " << s.getValueCart(0.,0.,x) << "  " << s.getErrorCart(0.,0.,x)<<endl;
    }
    ofstream dataOut("3d_legendre_data_rinv.dat");
    for (int i=0; i<62; ++i){
        double x = (double)i;
        dataOut << x << "  " << s.getValueSphr(0,0,x).real() << "  " << s.getErrorSphr(0,0,x).real()<<endl;
    }
    return true;
}

// ------------- harmonic_tester -------------
bool harmonic_tester(void){
    CCorrFtn3dSphr corr3d("K+", "K+", false, 2, "mini");
    parameterMap m;
    corr3d.Write(m);
    cout << m <<endl;
    corr3d.readTerms();
    parameterMap m0;
    corr3d[CSphericalHarmonicBasisFunction(0,0,true)].Write(m0);
    (corr3d.find(0,0,true))->second.Write(m0);
    corr3d.getItem(0,0,true).Write(m0);
    corr3d(0,0,true).Write(m0);
    cout << corr3d.getValueCart(10.,0.,0.) << " " << corr3d(0,0,true).getValue(10.)<<
    endl;
//    cout << m0 << endl;
    return true;
}

// ------------- object_io_tester -------------
bool object_io_tester(void){
    parameterMap com;
    ifstream f("sample_input/baselineModelData/testData_0_0_re.coral");
    f >> com;
    cout << com << endl;

    CCorrFtn3dSphr c;
    c.Read(com);
    
    for ( CCorrFtn3dSphr::iterator it=c.begin(); it!=c.end(); ++it) {
        cout  << it->first.l << " "  << it->first.m  << " " << boolalpha << it->first.realpart << endl;
    }     
    double r=10.;
    for (int l=0;l<=c.lmax;++l){
        for (int m=0;m<=l;++m){
            cout << l << ", " << m << ", "<<r<<": ";
            cout << c.getValueSphr(l,m,r);
            cout << " +/- ";
            cout << c.getErrorSphr(l,m,r) << endl;
        }
    }
    
    parameterMap tmp;
    c(0,0,true).Write(tmp);
    cout << tmp<< endl;
    
    return true;
}

// ------------- parameterMap_tester -------------
bool parameterMap_tester(void){
    parameterMap fake;
    parameter::set(fake,"fred",43);
    parameter::set(fake,"OK",string("fred"));
    string buf("insert\nlots of \ntext here!");
    parameter::set(fake,"big line",string("insert\nlots of \ntext here!"));
    parameter::set(fake,"a big line",buf);

    cout << fake << endl;

    string datafile("test_data/pmap/testData_0_0_re.coral");
    ifstream f(datafile.c_str());
    
    if (!f.good()) cerr << "\nOpening data file "<< datafile<<" failed!"<<endl;
    while (f.good()){
        parameterMap c;
        f >> c;
        cout << "\n"<<"---------------------------";
        cout << "\nParameters:\n";
        cout << c <<endl;
    }
    f.close();
    return true;
}
