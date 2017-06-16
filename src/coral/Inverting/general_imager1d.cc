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
#include "message.h"
#include "parametermap.h"
#include "general_imager1d.h"
#include "kernel_chooser.h"
#include "tnt_array1d_utils.h"
#include "linalg.h"

using namespace TNT;
using namespace std;

//--------------------------------------------
// Read/write to parameter map
//--------------------------------------------
//---------------- Read ----------------------------
bool CGeneralImager1d::Read( const parameterMap& m ){

//    parameterMap corrMap, souMap;
    ndata_corr = parameter::getI( m, "ndata_corr", ndata_corr );
    ndata_sou = parameter::getI( m,  "ndata_sou",  ndata_sou );
/*
    if ( m.find("source_settings")!=m.end() ){
        souMap  = parameter::getMap( m, "source_settings" );
        ndata_sou = parameter::getI(souMap,"ndata",ndata_sou);
        rmax = parameter::getD(souMap,"xmax",rmax);
        rmin = parameter::getD(souMap,"xmin",rmin);
    }
    if ( m.find("correlation_settings")!=m.end() ) {
        corrMap = parameter::getMap( m, "correlation_settings" );
        ndata_corr = parameter::getI(corrMap,"ndata",ndata_corr);
        qmax = parameter::getD(corrMap,"xmax",qmax);
        qmin = parameter::getD(corrMap,"xmin",qmin);
    }
*/
    rmax             = parameter::getD( m, "rmax", rmax );
    rmin             = parameter::getD( m, "rmin", rmin );
    qmax             = parameter::getD( m, "qmax", qmax );
    qmin             = parameter::getD( m, "qmin", qmin );
//    dq               = parameter::getD(corrMap,"dx",(qmax-qmin)/ndata_corr);
//    q0               = parameter::getD(corrMap,"xoffset",qmin);
//    bigQ             = parameter::getB(corrMap,"bigQ",false);
    if ( m.find("override_kernel") != m.end() ) {
        kernel_particle1 = parameter::getS( m, "kernel_particle1", kernel_particle1 );
        kernel_particle2 = parameter::getS( m, "kernel_particle2", kernel_particle2 );
    }
    return true;
}

//--------------- Write -----------------------------
// not used much except for diagnostics, so write everything!
bool CGeneralImager1d::Write( parameterMap& m ){

/*
    parameterMap corrMap, souMap;
    parameter::set(corrMap,"xmax",qmax);
    parameter::set(souMap,"xmax",rmax);
    parameter::set(corrMap,"xmin",qmin);
    parameter::set(souMap,"xmin",rmin);
    parameter::set(corrMap,"dx",dq);
    parameter::set(corrMap,"xoffset",q0);
    parameter::set(corrMap,"bigQ",bigQ);
    parameter::set(m,"source_settings",souMap);
    parameter::set(m,"correlation_settings",corrMap);
    parameter::set(corrMap,"ndata_corr",ndata_corr);
    parameter::set(souMap,"ndata_sou",ndata_sou);
*/
    parameter::set( m, "ndata_corr", ndata_corr );
    parameter::set( m, "ndata_sou",  ndata_sou );
    parameter::set( m, "rmax", rmax );
    parameter::set( m, "rmin", rmin );
    parameter::set( m, "qmax", qmax );
    parameter::set( m, "qmin", qmin );
    parameter::set( m, "override_kernel", true );
    parameter::set( m, "kernel_particle1", kernel_particle1 );
    parameter::set( m, "kernel_particle2", kernel_particle2 );
    return true;
}

//------------------- get_usable_data -------------------------
//! \brief Prepare the working copy of the correlation 
//! \param corrin The original correlation
//! This routine prepares the working copy of the correlation (CGeneratImager1d.corrwork)
//! corrwork has some important differences from the original:
//!    - bigQ is always false, so q = 0.5 (q1-q2)
//!    - the first data point (which is usually crap) is ignored
//!    - bad high q points are removed
//!    - 1 is subtracted off the l=0 terms in keeping with the convention of the Koonin-Pratt equation as implemented by Brown and Danielewicz
void CGeneralImager1d::get_usable_data(const CCorrFtn1dHisto& corrin){    

    // Determine which bins we keep
    int iqmin = std::max( corrin.whatBin( qmin ), 0 );
    int iqmax = std::min( corrin.whatBin( qmax ), corrin.ndata - 1 );
    ndata_corr = iqmax-iqmin + 1;
    if ( ndata_corr <= 1 ) {
        cout << "    No usable data in dataset, ndata_corr = "<<ndata_corr<<endl;
        return;
    }
    else if (ndata_corr == corrin.ndata) cout << "    Using all "          << ndata_corr << " are from this dataset"    << endl;
    else                                 cout << "    Using ndata_corr = " << ndata_corr << " data points from dataset" << endl;
    
    
    // Copy the initial correlation before messing with the copy
    corrwork.CopyState( corrin );
    
    // Fix the q to be q = 0.5*( q1 - q2 ), set the grid edges appropriately
    if ( corrwork.bigQ ){
        MESSAGE << "bigQ == true, resetting for corrwork" << ENDM_WARN;
        corrwork.bigQ = false;
        corrwork.xmin = corrin.leftBinEdge( iqmin ) / 2.0;
        corrwork.xmax = corrin.rightBinEdge( iqmax ) / 2.0;
    } else {
        corrwork.xmin = corrin.leftBinEdge( iqmin );
        corrwork.xmax = corrin.rightBinEdge( iqmax );
    }
    
    // Get copy of correlation, up to qmax, subtracting 1 when l=0.  
    Array1D< double > __data(   ndata_corr, 0.0 );
    Array1D< double > __uncert( ndata_corr, 0.0 );
    Array2D< double > __covmtx( ndata_corr, ndata_corr, 0.0 );   
    for ( int i=0; i<ndata_corr; ++i ){
        if ( l==0 ) __data[i] = corrin.data[i+iqmin] - 1.0; //corrin.data[ ndata_corr-1 ]; //- 1.0;
        else        __data[i] = corrin.data[i+iqmin];
        if ( corrin.covmtx_is_active ) {
            for ( int j=0; j<ndata_corr; ++j ) __covmtx[i][j] = corrin.covmtx[i+iqmin][j+iqmin];
            __uncert[i]=__covmtx[i][i];
        } else {
            __uncert[i] = corrin.uncert[i+iqmin];
            __covmtx[i][i] = corrin.uncert[i+iqmin] * corrin.uncert[i+iqmin];
        }
    }
    corrwork.ndata  = ndata_corr;
    corrwork.data   = __data;
    corrwork.uncert = __uncert;
    corrwork.covmtx = __covmtx;

    // 1st bin for l>0 always bogus, so we'll eliminate it if needed
    if ( ( l!=0 )&&( iqmin==0 ) ) { 
        for (int i=0;i<corrwork.ndata;++i){
            corrwork.covmtx[0][i]=0.0;
            corrwork.covmtx[i][0]=0.0;
        }
        corrwork.covmtx[0][0]=1.0;
        corrwork.data[0]=0.0;
        corrwork.uncert[0]=0.0;
    }
}

//------------------- set_kernel -------------------------
void CGeneralImager1d::set_kernel( const parameterMap& m ){

    kernelPtr = chooseKernel( kernel_particle1, kernel_particle2, m );
}

//------------------- set_kernel -------------------------
void CGeneralImager1d::set_kernel( const CKernel* _kernelPtr ){
    kernelPtr = _kernelPtr;
}

//------------------- initialize_source -------------------------
bool CGeneralImager1d::initialize_source( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m ){

    CObject1d* souPtr = dynamic_cast<CObject1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CObject1d"<<ENDM_FATAL;
    souPtr->CObject1d::CopyState( corrin );
    souPtr->Read(m);
    souout.particle1 = corrin.particle1;
    souout.particle2 = corrin.particle2;
    if ( ndata_corr <= 1 ) {
        set_no_data( souout );
        return false;
    }
    return true;
}

//------------------- initialize_correlation -------------------------
bool CGeneralImager1d::initialize_correlation( const CSourceFtnBase& souin, CCorrFtn1dHisto& corrout, const parameterMap& m ){

    const CObject1d* souPtr = dynamic_cast<const CObject1d*>(&souin);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CObject1d"<<ENDM_FATAL;
    corrout.CObject1d::CopyState( *souPtr );
    corrout.Read(m);
    corrout.particle1 = souin.particle1;
    corrout.particle2 = souin.particle2;
    return true;
}
