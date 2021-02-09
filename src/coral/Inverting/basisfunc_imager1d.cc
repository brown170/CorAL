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
#include "basisfunc_imager1d.h"

#include "message.h"
#include "constants.h"

#include "tnt_array1d.h"
#include "linalg.h"
#include "lsqrinvert.h"
#include "cubature_old.h"
#include "sf.h"
#include <cmath>

using namespace TNT;
using namespace std;

//#define __USE_GSL__

//#define __IMAGER_EPSILON 1e-10
#define __TYPICAL_SOURCE_SCALE__ 1e-4 //in fm^-3
//#define __TYPICAL_SOURCE_SCALE__ 1.0e-10 //in fm^-3
//#define __TYPICAL_SOURCE_SCALE__ 1.0e0 //in fm^-3

//---------------- Read ----------------------------
bool CBasisFuncImager1d::Read( const parameterMap& m ){

    MESSAGE << "CBasisFuncImager1d::Read"<<ENDM_INFO;
    
    CGeneralImager1d::Read( m );
    constrain_origin            = parameter::getB( m, "constrain_origin", false );
    constrain_rmin_zero         = parameter::getB( m, "constrain_rmin_zero", false );
    constrain_rmax_zero_slope   = parameter::getB( m, "constrain_rmax_zero_slope", false );
    constrain_rmax_zero         = parameter::getB( m, "constrain_rmax_zero", false );
    return true;
}

//--------------- Write -----------------------------
// not used much except for diagnostics, so write everything!
bool CBasisFuncImager1d::Write( parameterMap& m ){

    CGeneralImager1d::Write( m );
    parameter::set(m,"constrain_origin",constrain_origin);
    parameter::set(m,"constrain_rmin_zero",constrain_rmin_zero);
    parameter::set(m,"constrain_rmax_zero_slope",constrain_rmax_zero_slope);
    parameter::set(m,"constrain_rmax_zero",constrain_rmax_zero);
    return true;
}
    
//---------------- convertCorrelationToSource ----------------------------
//! function that manages the imaging
bool CBasisFuncImager1d::convertCorrelationToSource( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m, const CKernel* _kernelPtr ){

    CBasisFunctionExpansion1d* souPtr = dynamic_cast<CBasisFunctionExpansion1d*>(&souout);
    if (souPtr==NULL) 
        throw MESSAGE<<"Source argument must derive from CBasisFunctionExpansion1d"<<ENDM_FATAL;

    // Initialize this instance of the imager class
    cout <<"    Initializing imager, "<<flush;

    // ... Set internal parameters from corrin before reading map to establish defaults
    l                = corrin.l;
    ndata_corr       = corrin.ndata;
    qmax             = corrin.rightBinEdge(corrin.ndata-1);
    kernel_particle1 = corrin.particle1;
    kernel_particle2 = corrin.particle2;

    // ... Read rest of parameters from map (override as needed)
    Read( parameter::getMap( m, "imager_settings" ) );

    // ... Override parameters with data from correlation 
    qmax             = std::min(corrin.rightBinEdge(corrin.ndata-1),qmax);
    qmin             = std::max(corrin.leftBinEdge(0),qmin);
    dq               = corrin.binWidth(0);
    q0               = qmin+dq/2.;

    // Set the kernel
    cout <<"kernel, "<<flush;
    if ( _kernelPtr == NULL ) set_kernel( parameter::getMap( m, "kernel_settings" ) );
    else  set_kernel( _kernelPtr );

    // Get the usable data.
    cout <<"and working copy of data."<<endl;
    get_usable_data( corrin );

    // Initialize the output source.     
    if ( !initialize_source( corrwork, souout, parameter::getMap(m,"source_settings") ) ) return false;
    ndata_sou = souPtr->ndata; //min( , ndata_corr )
    cout << "    Source will have "<<ndata_sou<<" coefficients"<<endl;
    
    // Some simple checks on the correlation and the source 
    if ( ndata_corr <= 1 ) throw MESSAGE<<"Correlation has no data, ndata_corr = "<<ndata_corr<<ENDM_FATAL;
    if ( souPtr->ndata <= 1 ) throw MESSAGE<<"Source has no space for data, ndata_sou = "<<ndata_sou<<ENDM_FATAL;

    // Now do the imaging
    if ( kmtx != CKernelMatrix( l, ndata_corr, q0, dq, ndata_sou, true, kernel_particle1, kernel_particle2 ) ) 
    {
        cout << "    Generating kmtx:"<<endl;
        set_kmtx( corrwork, *souPtr, m );
    }  
        
    imageit( *souPtr );
    return true;
}

//------------------ convertSourceToCorrelation --------------------------
//! function that manages the unimaging
bool CBasisFuncImager1d::convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn1dHisto& corrout, const parameterMap& m, const CKernel* _kernelPtr ){

    const CBasisFunctionExpansion1d* souPtr = dynamic_cast<const CBasisFunctionExpansion1d*>(&souin);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisFunctionExpansion1d"<<ENDM_FATAL;
    
    // Set internal parameters from souin before reading map to establish defaults
    l                = souPtr->l;
    ndata_sou        = souPtr->ndata;
    rmax             = souPtr->getRightSupport( souPtr->ndata - 1 );
    rmin             = souPtr->getLeftSupport( 0 );
    kernel_particle1 = souin.particle1;
    kernel_particle2 = souin.particle2;

    // Read rest of parameters from map (override as needed)
    Read( parameter::getMap( m, "imager_settings" ) );

    // Set the kernel
    if ( _kernelPtr == NULL ) set_kernel( parameter::getMap( m, "kernel_settings" ) );
    else  set_kernel( _kernelPtr );

    // Initialize the output correlation    
    if ( !initialize_correlation( souin, corrout, parameter::getMap(m,"correlation_settings") ) ) return false;
    if ( ndata_sou <= 1 ) throw MESSAGE<<"Source has no data, ndata_sou = "<<ndata_sou<<ENDM_FATAL;
    if ( ndata_corr <= 1 ) ndata_corr = max( corrout.ndata, ndata_sou );
    corrout.covmtx_is_active=true;
    corrout.setDim( ndata_corr );
    corrout.setFixedWidthBins( dq, q0 );

    // Unimage
    if ( kmtx != CKernelMatrix( l, ndata_corr, q0, dq, ndata_sou, true, kernel_particle1, kernel_particle2 ) ) {
       cout << "    Generating kmtx:"<<endl;
       set_kmtx( corrout, *souPtr, m );
    }
    unimageit( *souPtr, corrout );
    return true;
}

//--------------------- unimageit ---------------------------
void CBasisFuncImager1d::unimageit(const CBasisFunctionExpansion1d& souin, CCorrFtn1dHisto& corrout){

    // Compute the correlation
    corrout.data=( 1.0 / __TYPICAL_SOURCE_SCALE__ ) * kmtx * souin.data;

    // Figure out the source covariance matrix
    Array2D< double > covmtx;
    if ( souin.covmtx_is_active ) covmtx = souin.covmtx;
    else covmtx = makeDiagSquared( souin.uncert );
    
    // Compute the correlation covariance matrix
    corrout.covmtx=( 1. / ( __TYPICAL_SOURCE_SCALE__*__TYPICAL_SOURCE_SCALE__ ) ) * matmult( matmult( kmtx, covmtx ), transpose( kmtx ) );
    corrout.syncUncert();
    
    // For l=0, have to add back 1.0
    if ( l==0 ) {
        for ( int i=0; i<corrout.data.dim(); i++ ) corrout.data[i] = corrout.data[i] + 1.0;
    }
}

//--------------------- imageit ---------------------------
//! Code to image a 1d source, with equality constraints
double CBasisFuncImager1d::imageit(CBasisFunctionExpansion1d& souout){

    // Set up the constraints
    set_constraints(souout);

cout << "kmtx:"<<kmtx<<endl;
cout <<"conmtx:"<<conmtx<<endl;
cout<< "convec"<<convec<<endl;
//parameterMap junk;
//corrwork.Write( junk );
//souout.Write( junk );
//cout << "corrwork:"<<junk << endl;

    // Now do the inversion
    cout << "  Inverting..."<<endl;
#undef CHEAT
#ifdef CHEAT
    souout.data = svdinvert( kmtx ) * corrwork.data;
    souout.uncert = svdinvert( kmtx ) * corrwork.uncert;
    souout.covmtx_is_active = false;
#else
    if (num_constraints>0) {
        cout << "    Using constraints"<<endl;
        CLSqrInvertSVDLagrange inverter( kmtx, conmtx, convec );
        inverter.solve( corrwork.data, corrwork.covmtx );
        souout.data   = inverter.model();
        souout.covmtx = inverter.covmodel();
        cout << "    Lagrange Multiplier vector:\n";
        cout << inverter.lagrange_multipliers() << endl;
    } else 
        LeastSquaresInvert( corrwork.data, corrwork.covmtx, kmtx, souout.data, souout.covmtx );
#endif

    // rescale the source back (we rescaled the kmtx to get O(1))
    souout.data   = __TYPICAL_SOURCE_SCALE__*souout.data;
    souout.covmtx = __TYPICAL_SOURCE_SCALE__*__TYPICAL_SOURCE_SCALE__*souout.covmtx;
    souout.covmtx_is_active = true;

    // Sync up the uncertainty estimates, they're the sqrt(diag(covmtx))
    cout << "    Synching uncertainty with covariance" << endl;
    souout.syncUncert();
  
//cout << "S(0): " <<souout.getValue( 0.0 ) << "+/-" << souout.getError( 0.0 ) << endl;
    
    // Compute final chi^2 (obvious thing to do)
    double chi2 = 0.0, tmp = 0.0;
    Array1D<double> residual = ( 1.0/__TYPICAL_SOURCE_SCALE__ ) * kmtx * souout.data - corrwork.data;
    for ( int i=0; i<corrwork.ndata; ++i ){
        tmp = residual[i] / sqrt( corrwork.covmtx[i][i] );
        chi2 += tmp * tmp;
    }
    
    if (isnan(chi2)) throw MESSAGE << "Imaging failed, chi^2 is NaN!"<< ENDM_FATAL;
    if (chi2/corrwork.ndata>1e3) MESSAGE<<"Terrible chi^2, something is probably wrong!"<<ENDM_SEVERE;
    if (chi2/corrwork.ndata<1e-3) MESSAGE<<"Way too good chi^2, you're over resolving something!"<<ENDM_SEVERE;
    cout << "    Final chi^2 = " << chi2 << endl;
    cout << "    Final chi^2/NDF = " << chi2/corrwork.ndata << endl;
    
    return chi2;
}

//------------------- set_no_data -------------------------
bool CBasisFuncImager1d::set_no_data( CSourceFtnBase& souout ){

    cout << "    Zeroed data, skipping the imaging...\n";
    Array1D<double> tmpcoef(ndata_sou,0.0);
    Array2D<double> tmpcov(ndata_sou,ndata_sou,0.0);
    CBasisFunctionExpansion1d* souPtr = dynamic_cast< CBasisFunctionExpansion1d* >( &souout );
    if (souPtr==NULL) throw MESSAGE << "Source argument must derive from CBasisFunctionExpansion1d" << ENDM_FATAL;
    souPtr->data   = tmpcoef;
    souPtr->covmtx = tmpcov;
    souPtr->ndata  = ndata_sou;
    souPtr->xmin = 0.0;
    souPtr->xmax = rmax;
    souPtr->syncUncert();
    souPtr = NULL;
    return true;
}

//------------------- initialize_source -------------------------
bool CBasisFuncImager1d::initialize_source( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m ){

    CBasisFunctionExpansion1d* souPtr = dynamic_cast<CBasisFunctionExpansion1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisFunctionExpansion1d"<<ENDM_FATAL;
    return CGeneralImager1d::initialize_source(corrin,souout,m);
}

//----------------- set_constraints ---------------------------
void CBasisFuncImager1d::set_constraints( const CBasisFunctionExpansion1d& souout ){

    if ( ndata_sou != souout.ndata ) 
        throw MESSAGE << "ndata_sou ("<<ndata_sou<<") != souout.ndata ("<<souout.ndata<<")!!!" <<ENDM_FATAL;

    cout << "  Setting up constraints"<<endl;

    // count constraints
    num_constraints=0;
    if (constrain_origin) {
        if (souout.l == 0) {num_constraints+=1;} else {num_constraints+=2;}
    }
    if (constrain_rmax_zero_slope) num_constraints+=1;
    if (constrain_rmax_zero) num_constraints+=1;
    if (num_constraints>0) 
        cout<<"    Initializing "<<num_constraints<<" valid constraints\n";
    else {
        cout<<"    No constraints to initialize\n";
        return;
    }

    // allocate constraint matrix and vector
    Array2D<double> _conmtx(num_constraints,souout.ndata);
    Array1D<double> _convec(num_constraints,0.0);
    conmtx=_conmtx;
    convec=_convec;   
    int icons=0; // constraint index counter

    // build r=0 constraints
    // dS_lm(0)/dr = 0 *always*
    if (constrain_origin) {
        cout << "      Constraining dS_lm(0)/dr = 0"<<endl;
        for (int i=0;i<souout.ndata;++i){conmtx[icons][i]=souout.basisFunction(0.0,i,1);}
        icons++;
    }

    // S_lm(0) = 0 for all lm!=00
    if ((souout.l!=0)&&constrain_origin) {
        cout << "      Constraining S_lm(0) = 0 for l > 0"<<endl;
        for (int i=0;i<souout.ndata;++i){conmtx[icons][i]=souout.basisFunction(0.0,i,0);}
        icons++;
    }

    // set S_lm(rmin) = 0 
    if ((souout.l!=0)&&constrain_rmin_zero) {
        cout << "      Constraining S_lm(rmin) = 0"<<endl;
        for (int i=0;i<souout.ndata;++i){conmtx[icons][i]=souout.basisFunction(rmin,i,0);}
        icons++;
    }

    // build S_lm(rmax) = 0 constraint
    if (constrain_rmax_zero) {
        cout << "      Constraining S_lm(rmax) = 0"<<endl;
        for (int i=0;i<souout.ndata;++i){conmtx[icons][i]=souout.basisFunction(rmax,i,0);}
        icons++;
    }

    // build dS_lm(rmax)/dr = 0 constraint
    if (constrain_rmax_zero_slope) {
        cout << "      Constraining dS_lm(rmax)/dr = 0"<<endl;
        for (int i=0;i<souout.ndata;++i){conmtx[icons][i]=souout.basisFunction(rmax,i,1);}
        icons++;
    }
}

//----------------- set_kmtx ---------------------------
void CBasisFuncImager1d::set_kmtx( const CCorrFtn1dHisto& corrin, const CBasisFunctionExpansion1d& souout, const parameterMap& m ){

    sourcePtr = &souout;
    CKernelMatrix kmtx_tmp(corrin.l, ndata_corr, q0, dq, ndata_sou, true, kernel_particle1, kernel_particle2);

    CIntegrateCubature integrate;
    double dq, qmin, qmax, rmin, rmax;
    integrate.set_ndim(2);
    integrate.set_fdim(1);
    __l=corrin.l;
    for (int i=0;i<kmtx_tmp.dim1();++i){
        dq   = corrin.binWidth(i);
        qmin = corrin.getLeftSupport(i);
        qmax = corrin.getRightSupport(i);
        integrate.set_limit( 0, qmin, qmax );
        for ( int j=0; j<kmtx_tmp.dim2(); ++j ) {
            rmin = sourcePtr->getLeftSupport(j);
            rmax = sourcePtr->getRightSupport(j);
            __j  = j;
            integrate.set_limit( 1, rmin, rmax );
            if ( qmax <= qmin ) kmtx_tmp[i][j] = 0.0;
            else {
                integrate.compute( kp_integrand, this );
                kmtx_tmp[i][j] = __TYPICAL_SOURCE_SCALE__ * integrate.value / dq;
            }
        }
    }

    kmtx      = kmtx_tmp;
    sourcePtr = NULL;
    cout <<"    The kmtx has been generated!"<<endl;
}

//----------------- kp_integrand ---------------------------
int CBasisFuncImager1d::kp_integrand( unsigned int ndim, const double* x, void* classptr, unsigned fdim, double *fval ){//changed return type of function from double to int

    CBasisFuncImager1d *cls = (CBasisFuncImager1d*)classptr;
//    return 1.0;

    return 4.0 * PI * x[1] * x[1] *
//Bessel::jn(cls->__l,2.0*x[0]*x[1]/197.3269602)* // for debugging
        ( cls->kernelPtr->GetValue( cls->__l, x[0], x[1] ) )*
        ( cls->sourcePtr->basisFunction( x[1], cls->__j, 0 ) );
};

