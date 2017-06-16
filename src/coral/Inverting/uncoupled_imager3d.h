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
#ifndef UNCOUPLED_IMAGER3D_H
#define UNCOUPLED_IMAGER3D_H

#include <string>
#include "parametermap.h"
#include "kernel.h"
#include "corr3d_ylm.h"
#include "sou3d_ylm.h"
#include "soubase.h"
#include "sou1d_bsplines.h"
#include "sou1d_legendre.h"
#include "sou1d_laguerre.h"
#include "sou1d_hermite.h"
#include "sou1d_chebyshev.h"
#include "sou1d_histo.h"
#include "general_imager3d.h"
#include "bspline_imager1d.h"
#include "basisfunc_imager1d.h"

using namespace std;
using namespace TNT;

//------------------------------------------
// Main imaging code
//------------------------------------------
template< class TSource1d, class TImager1d >
class CUncoupledImager3d: public CGeneralImager3d {

public:
    
    // ----------- CUncoupledImager3d -----------------
    //! Constructor
    CUncoupledImager3d( void ): CGeneralImager3d(){}

    // ----------- ~CUncoupledImager3d -----------------
    //! Destructor
    ~CUncoupledImager3d( void ){}

    // ----------- Read -----------------
    //! Read from parameter map
    bool Read( const parameterMap& m ){return CGeneralImager3d::Read(m);}

    // ----------- Write -----------------
    //! write to parameter map
    bool Write( parameterMap& m ){return CGeneralImager3d::Write(m);}

    // ----------- convertCorrelationToSource -----------------
    //! function that actually manage the imaging
    bool convertCorrelationToSource( const CCorrFtn3dSphr& corrin, CSourceFtnBase& souout, const parameterMap& m ){

        CSourceFtn3dSphr< TSource1d >* souPtr = dynamic_cast< CSourceFtn3dSphr< TSource1d >* >(&souout);
        if (souPtr == NULL) throw MESSAGE << "Source argument must derive from CSourceFtn3dSphr<TSource1d>"<<ENDM_FATAL;

        parameterMap sou_m = parameter::getMap( m, "source_settings" );
        souPtr->lmax = parameter::getI( sou_m, "lmax", corrin.lmax );
        souPtr->particle1 = parameter::getS( sou_m, "particle1", corrin.particle1 );
        souPtr->particle2 = parameter::getS( sou_m, "particle2", corrin.particle2 );
        souPtr->storage_directory = parameter::getS( sou_m, "storage_directory", "." );
        souPtr->skip_odd_l = parameter::getB( sou_m, "skip_odd_l", souout.particle1 == souout.particle2 );

        // Set the kernel
        if ( kernelPtr == NULL ) set_kernel( parameter::getMap( m, "kernel_settings" ) );

        bool result=true;
        TImager1d imager;
        for ( CCorrFtn3dSphr::const_iterator it=corrin.begin(); it!=corrin.end(); ++it ){
            TSource1d souterm;
            cout << "Imaging .... term: " << it->first.termName() << endl;
            imager.convertCorrelationToSource( it->second, souterm, m, kernelPtr );
            souPtr->insert( make_pair( it->first, souterm ) );
        }
        return result;
    }
    
    // ----------- convertSourceToCorrelation -----------------
    //! function that actually manage the unimaging
    bool convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn3dSphr& corrout, const parameterMap& m ){
        const CSourceFtn3dSphr<TSource1d>* souPtr = dynamic_cast< const CSourceFtn3dSphr<TSource1d>* >(&souin);
        if (souPtr == NULL) throw MESSAGE << "Source argument must derive from CSourceFtn3dSphr<TSource1d>"<<ENDM_FATAL;

        parameterMap cor_m = parameter::getMap( m, "correlation_settings" );
        corrout.lmax=parameter::getI(cor_m,"lmax",souPtr->lmax);
        corrout.bigQ=parameter::getB(cor_m,"bigQ",false);
        corrout.particle1=parameter::getS(cor_m,"particle1",souin.particle1);
        corrout.particle2=parameter::getS(cor_m,"particle2",souin.particle2);
        corrout.storage_directory=parameter::getS(cor_m,"storage_directory",".");
        corrout.skip_odd_l=parameter::getB(cor_m,"skip_odd_l",corrout.particle1==corrout.particle2);

        // Set the kernel
        if ( kernelPtr == NULL ) set_kernel( parameter::getMap( m, "kernel_settings" ) );

        bool result=true;
        TImager1d imager;
        for ( 
            typename CSourceFtn3dSphr<TSource1d>::const_iterator it=souPtr->CSourceFtn3dSphr<TSource1d>::begin(); // quite a mouthful!
            it!=souPtr->CSourceFtn3dSphr<TSource1d>::end(); 
            ++it 
        ){
            cout <<"Unimaging .... term: "<<it->first.termName()<<endl;
            CCorrFtn1dHisto corrterm;
            imager.convertSourceToCorrelation( it->second, corrterm, m, kernelPtr );
            corrout.insert( make_pair( it->first, corrterm ) );
        }
        return result;
    }
        
    // ----------- get_usable_data -----------------
    //void get_usable_data( const CCorrFtn3dSphr& corrin );

    // ----------- initialize_source -----------------
    //! unused
    bool initialize_source( const CCorrFtn3dSphr& corrin, CSourceFtnBase& souout ){ return true; }
    
    // ----------- initialize_correlation -----------------
    //! unused
    bool initialize_correlation( const CSourceFtnBase& souin, CCorrFtn3dSphr& corrout ){ return true; }
    
    // ----------- set_no_data -----------------
    //! unused
    bool set_no_data( CSourceFtnBase& souout ){ return true; }
};

typedef CUncoupledImager3d< CSourceFtn1dBSpline,       CBasisSplineImager1d > UncoupledBasisSplineImager3d;
typedef CUncoupledImager3d< CSourceFtn1dHisto,         CBasisSplineImager1d > UncoupledHistoImager3d;
typedef CUncoupledImager3d< CSourceFtn1dLegendrePoly,  CBasisFuncImager1d >   UncoupledLegendrePolyImager3d;
typedef CUncoupledImager3d< CSourceFtn1dLaguerrePoly,  CBasisFuncImager1d >   UncoupledLaguerrePolyImager3d;
typedef CUncoupledImager3d< CSourceFtn1dHermitePoly,   CBasisFuncImager1d >   UncoupledHermitePolyImager3d;
typedef CUncoupledImager3d< CSourceFtn1dChebyshevPoly, CBasisFuncImager1d >   UncoupledChebyshevPolyImager3d;

#endif
