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
#include "kernel.h"
#include "crombergintegrator.h"
#include "constants.h"

using namespace std;

//-------------------------------------------------------------------------
// Function prototypes
//-------------------------------------------------------------------------
void plotSource( yasper::ptr< CBasisFunctionExpansion1d > pExpansion, const parameterMap& inMap );
void plotCorrelation( yasper::ptr< CBasisFunctionExpansion1d > pExpansion, CKernel* pKernel, const parameterMap& inMap );
void getHelp(void);
int main(int argc, char* argv[]);

//-------------------------------------------------------------------------
// Classes
//-------------------------------------------------------------------------
class CBFPLOTIntegrand{
  public:
    CBFPLOTIntegrand( yasper::ptr< CBasisFunctionExpansion1d > _pExpansion, CKernel* _pKernel, int _i, int _l, double _q ):
        pExpansion(_pExpansion), pKernel(_pKernel), index(_i), l(_l), q(_q){}
    yasper::ptr< CBasisFunctionExpansion1d > pExpansion;
    CKernel* pKernel;
    int index, l;
    double q;
    static double f( void* instance, double r ){
        CBFPLOTIntegrand* inst = static_cast< CBFPLOTIntegrand* >(instance);
        return 4.0 * PI * r * r *
               inst->pKernel->GetValue( inst->l, inst->q, r ) *
               inst->pExpansion->basisFunction( r, inst->index );
    }
};

//-------------------------------------------------------------------------
// Function implementations
//-------------------------------------------------------------------------

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"\n";
    cout<<"\n";
    cout<<"\nUsage: bfplot <option> [inputFile.dat]"<<endl;
    cout<<"    options -h, -help, --help all print this message, then quit"<<endl;
    cout<<"    -legendre       Print out Legendre polynomial basis" <<endl;
    cout<<"    -bspline        Print out Basis Spline basis"   <<endl;
    cout<<"    -laguerre       Print out Laguerre function basis"  <<endl;
    cout<<"    -chebyshev      Print out Chebyshev polynomial basis" <<endl;
    cout<<"    -histogram      Print out Histogram basis" <<endl;
    cout<<"    -hermite        Print out Hermite function basis"   <<endl;
    cout<<"    -convolve       Convolve basis funcs with kernel too"   <<endl;
    exit(0);
}

// ------------- main -------------
int main(int argc, char* argv[]){

    cout << "*** Widget to PLOT the Basis Functions (BFPLOT) using CorAL ***"<<endl;
    cout << endl;

    bool got_file = false;
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
        if (sarg=="-convolve") convolve_with_kernel=true;
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

    yasper::ptr<CBasisFunctionExpansion1d> pExpansion;
    for ( vector<string>::iterator it=modeList.begin(); it!=modeList.end(); ++it)
    {
        if      ( *it == "-legendre" )  { pExpansion = new CLegendrePolynomialExpansion1d; }
        else if ( *it == "-bspline" )   { pExpansion = new CBasisSpline1d; }
        else if ( *it == "-laguerre" )  { pExpansion = new CLaguerrePolynomialExpansion1d; }
        else if ( *it == "-chebyshev" ) { pExpansion = new CChebyshevPolynomialExpansion1d; }
        else if ( *it == "-histogram" ) { pExpansion = new CHistogram1d; }
        else if ( *it == "-hermite" )   { pExpansion = new CHermiteFunctionExpansion1d; }
        else  MESSAGE<<"Unknown mode: '"<<*it<<ENDM_WARN;
    }
    pExpansion->Read( inMap );
    plotSource( pExpansion, inMap );

    if (inMap.hasKey("particle1")&&inMap.hasKey("particle2")&&inMap.hasKey("param_filename")&&convolve_with_kernel) {
        string p1 = parameter::getS(inMap, "particle1", "pi+");
        string p2 = parameter::getS(inMap, "particle2", "pi+");
        CKernel* pKernel = chooseKernel( p1, p2, inMap );
        plotCorrelation( pExpansion, pKernel, inMap );
    }

    return 0;
}

// ------------- plotSource -------------
void plotSource( yasper::ptr< CBasisFunctionExpansion1d > pExpansion, const parameterMap& inMap ){
    string outFile = parameter::getS(inMap,"basis_func_plot_file","output_r_basis_funcs.dat");
    cout << "    Plotting basis functions in r-space to " << outFile << endl;
    ofstream fout(outFile.c_str());
    for (int i=0; i<pExpansion->ndata; ++i) {
        fout << endl;
        for (int n=0;n<100;++n) {
            double r = (pExpansion->getLeftSupport(i)-pExpansion->getRightSupport(i))*double(n)/100.0 + pExpansion->getRightSupport(i);
            fout << r << " " << pExpansion->basisFunction(r,i)<<endl;
        }
    }

}

// ------------- plotCorrelation -------------
void plotCorrelation( yasper::ptr< CBasisFunctionExpansion1d > pExpansion, CKernel* pKernel, const parameterMap& inMap  ){
    string outFile = parameter::getS(inMap,"convoluted_basis_func_plot_file","output_q_basis_funcs.dat");
    cout << "    Plotting basis functions in q-space to " << outFile << endl;
    ofstream fout(outFile.c_str());
    double qmin=parameter::getD(inMap,"qmin",0.0);
    double qmax=parameter::getD(inMap,"qmax",100.0);
    CBFPLOTIntegrand integrand( pExpansion, pKernel, 0, parameter::getI(inMap,"l",0), 0.0 );
    CRombergIntegrator integrator( static_cast<void*>(&integrand), CBFPLOTIntegrand::f );
    for (int i=0; i<pExpansion->ndata; ++i) {
        fout << endl;
        integrand.index = i;
        for (int n=0;n<100;++n) {
            double q = (qmax-qmin)*double(n)/100.0 + qmin;
            integrand.q = q;
            fout << q << " " << integrator.compute(pExpansion->getLeftSupport(i),pExpansion->getRightSupport(i))<<endl;
        }
    }

}
