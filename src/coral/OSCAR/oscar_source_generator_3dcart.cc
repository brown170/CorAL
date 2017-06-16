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
#include "oscar_source_generator_3dcart.h"
#include "constants.h"
#include "tnt_array1d_utils.h"
#include "misc.h"
#include "cheezyparser.h"



// ------------------ Read ------------------
bool COSCARSourceGenerator3dCart::Read( const parameterMap& m ){
    bool result = COSCARSourceGeneratorBase::Read(m);
    fold[iSide] = parameter::getB( m, "fold_side", fold[iSide] );
    fold[iOut]  = parameter::getB( m, "fold_out",  fold[iOut]  );
    fold[iLong] = parameter::getB( m, "fold_long", fold[iLong] );
    return result;
}


// ------------------ Write ------------------
bool COSCARSourceGenerator3dCart::Write( parameterMap& m){
    bool result = COSCARSourceGeneratorBase::Write(m);
    parameter::set( m, "fold_side", fold[iSide] );
    parameter::set( m, "fold_out",  fold[iOut]  );
    parameter::set( m, "fold_long", fold[iLong] );
    return result;
}


// ------------------ addOnePair ------------------
bool COSCARSourceGenerator3dCart::addOnePair( const COSCARLine& p1, const COSCARLine& p2 ){
    // bin up 3D source
    double r[3];
    getSideOutLong( r_cm, r );
    // Generate the unnormalized source histogram
    for ( int i=0; i<3; ++i ) if ( fold[i] ) r[i] = abs( r[i] );
    // Get some stats on the r_cm's
    maxS = std::max( r[iSide], maxS );
    minS = std::min( r[iSide], minS );
    maxO = std::max( r[iOut],  maxO );
    minO = std::min( r[iOut],  minO );
    maxL = std::max( r[iLong], maxL );
    minL = std::min( r[iLong], minL );
    int iBin = result.whatBin( r[0], r[1], r[2] );
    if ( iBin >= 0 && iBin < result.ndata ) {
        result.data[iBin]   += 1;
        result.uncert[iBin] += 1;
    }
    return true;
}


// ------------------ postProcessPairs ------------------
bool COSCARSourceGenerator3dCart::postProcessPairs( void ){
    cout << "    Normalizing 3D source function...\n";
    double integral=0.0, dintegral=0.0;
    double flip_factor=1.0;
    for (int i=0;i<3;++i) if (fold[i]) flip_factor*=2.0;
    double rweight=result.binVolume()*flip_factor; ///;
    result.covmtx_is_active=false;
    for (int i=0;i<result.ndata;++i) {
        result.data[i]=result.data[i]/double(totalPairs)/rweight;
        integral += result.data[i]*rweight;
    }
    for (int i=0;i<result.ndata;++i) {
        double u = result.uncert[i]/double(totalPairs)/rweight/rweight-result.data[i]*result.data[i];
        if ( u < -1e-6 ) throw MESSAGE <<"Bad u: "<<u<<ENDM_FATAL;
        result.uncert[i]=sqrt(abs(u/double(totalPairs)));
        dintegral += pow(result.uncert[i]*rweight,2.0);
    }
    cout << "    Some parameters from the legal pairs:\n";
    cout << "        maxS:" << maxS << ", minS:" << minS << "\n";
    cout << "        maxO:" << maxO << ", minO:" << minO << "\n";
    cout << "        maxL:" << maxL << ", minL:" << minL << "\n";
    cout << "        Number of good pairs: " << totalPairs << endl;
    cout << "    Cross-checking integral of Cartesian binned source function..."<<endl;
    cout << "        Integral of source: " << integral << "+/-" << sqrt(dintegral) << endl;
    cout << "        'Flip factor' due to assumed symmetries: " << flip_factor << endl;
    return true;
}

// ------------------ generateSource ------------------
CSourceFtn3dHisto COSCARSourceGenerator3dCart::generateSource( vector<COSCARLine> plist, const parameterMap& m ){

    // Initialize the source we are generating
    parameterMap souMap; 
    souMap = parameter::getMap( m,"source_settings" );
    result.Read(souMap);

    // Initialize this instance of the source generator
    Read(m);
    
    // Filter the pairs
    filterParticles( plist );
    
    // Now build the source
    accumulatePairs();
    
    return result;
}
