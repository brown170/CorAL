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
#include "oscar_correlation_generator_3dcart.h"
#include "constants.h"
#include "tnt_array1d_utils.h"
#include "misc.h"
#include "cheezyparser.h"



// ------------------ Read ------------------
bool COSCARCorrelationGenerator3dCart::Read( const parameterMap& m ){
    bool result = COSCARCorrelationGeneratorBase::Read(m);
    fold[iSide] = parameter::getB( m, "fold_side", fold[iSide] );
    fold[iOut]  = parameter::getB( m, "fold_out",  fold[iOut]  );
    fold[iLong] = parameter::getB( m, "fold_long", fold[iLong] );
    return result;
}


// ------------------ Write ------------------
bool COSCARCorrelationGenerator3dCart::Write( parameterMap& m){
    bool result = COSCARCorrelationGeneratorBase::Write(m);
    parameter::set( m, "fold_side", fold[iSide] );
    parameter::set( m, "fold_out",  fold[iOut]  );
    parameter::set( m, "fold_long", fold[iLong] );
    return result;
}


// ------------------ addOnePair ------------------
bool COSCARCorrelationGenerator3dCart::addOnePair( const COSCARLine& p1, const COSCARLine& p2 ){
    // This is a good pair, let's bin it up the source function(s)...
    // bin up 3D source
    double r[3];
    double weight;
    double q_dot_r; 
    double costheta; 
    getSideOutLong( r_cm, r );
    
    if (accumulation_mode=="default"){
        double q[3];
        getSideOutLong( q_cm, q );
        q_dot_r  = (q[0]*r[0]+q[1]*r[1]+q[2]*r[2]);
        costheta = q_dot_r/qinv/rinv;
        // Get some stats on the r_cm's
        maxS = std::max( q[iSide], maxS );
        minS = std::min( q[iSide], minS );
        maxO = std::max( q[iOut],  maxO );
        minO = std::min( q[iOut],  minO );
        maxL = std::max( q[iLong], maxL );
        minL = std::min( q[iLong], minL );
        // Generate the unnormalized source histogram
        for (int i=0;i<3;++i) if (fold[i]) q[i]=abs(q[i]);
        int iBin = result.whatBin(q[0],q[1],q[2]);
        if (iBin>=0 && iBin<result.ndata) {
            weight = (kernelPtr->GetPsiSquared(qinv, rinv, costheta));
            result.data[iBin]   += weight;
            result.uncert[iBin] += weight*weight;
            pairCount[iBin]     += 1;
        }

    } else {
        int iBin;
        for (int i = 0; i < 10; ++i ){
            iBin = static_cast<int>( ranGen.iran(result.ndata) );
            Array1D< double > q = result.binCenter( iBin );
            q_dot_r = (q[0]*r[0]+q[1]*r[1]+q[2]*r[2]);
            double new_qinv=sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]);
            costheta = q_dot_r/rinv/new_qinv;
            maxS=std::max(q[0],maxS);
            minS=std::min(q[0],minS);
            maxO=std::max(q[1],maxO);
            minO=std::min(q[1],minO);
            maxL=std::max(q[2],maxL);
            minL=std::min(q[2],minL);
            weight = (kernelPtr->GetPsiSquared(new_qinv, rinv, costheta));
            result.data[iBin]   += weight;
            result.uncert[iBin] += weight*weight;
            pairCount[iBin]     += 1;
        }
    }
    return true;
}


// ------------------ postProcessPairs ------------------
bool COSCARCorrelationGenerator3dCart::postProcessPairs( void ){
    cout << "    Normalizing 3D correlation function...\n";
    result.covmtx_is_active=false;
    for (int i=0;i<result.ndata;++i) {
        if ( pairCount[i] > 1 ) {
            result.data[i]   = 1.0+result.data[i]/double(pairCount[i]);
            result.uncert[i] = sqrt(abs(result.uncert[i]/double(pairCount[i]) - result.data[i]*result.data[i])/double(pairCount[i]));
        } else {
            result.data[i]   = 1.0;
            result.uncert[i] = 1.0;
        }
    }
    cout << "    Some parameters from the legal pairs:\n";
    cout << "        maxS:" << maxS << ", minS:" << minS << "\n";
    cout << "        maxO:" << maxO << ", minO:" << minO << "\n";
    cout << "        maxL:" << maxL << ", minL:" << minL << "\n";
    return true;
}

// ------------------ generateCorrelation ------------------
CCorrFtn3dHisto COSCARCorrelationGenerator3dCart::generateCorrelation( vector<COSCARLine> plist, const parameterMap& m ){

    // Initialize this instance of the correlation generator
    Read(m);
    
    // Initialize the correlation we are generating
    parameterMap corrMap; 
    corrMap = parameter::getMap( m,"correlation_settings" );
    result.Read(corrMap);

    // Get the controls for the kernel we'll use to generate the correlation from the source
    parameterMap kernelMap;
    kernelMap = parameter::getMap(m, "kernel_settings");
    set_kernel(kernelMap);
 
    // Filter the pairs
    filterParticles( plist );
    
    // Setup the pair counter
    TNT::Array1D<int> newPairCount(result.ndata,0);
    pairCount = newPairCount;
    
    // Now build the source
    accumulatePairs();
    
    return result;
}
