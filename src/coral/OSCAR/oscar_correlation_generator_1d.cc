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
#include "oscar_correlation_generator_1d.h"
#include "constants.h"
#include "tnt_array1d_utils.h"
#include "misc.h"
#include "message.h"
#include "sf.h"


// ------------------ addOnePair ------------------
bool COSCARCorrelationGenerator1d::addOnePair( const COSCARLine& p1, const COSCARLine& p2 ){
    double r[3];
    double weight;
    getSideOutLong( r_cm, r );
    if (accumulation_mode=="default"){
        int iBin = result.whatBin(qinv,false);
        if (iBin>=0 && iBin<result.ndata) {
            if ( result.sphrHarm ) {
                if (result.realpart) weight =   SpherHarmonics::ReYlm(result.l,result.m,r[0],r[1],r[2]);
                else                 weight = - SpherHarmonics::ImYlm(result.l,result.m,r[0],r[1],r[2]);
                weight *= (kernelPtr->GetValue(result.l, qinv, rinv));
                result.data[iBin]   += weight;
                result.uncert[iBin] += weight*weight;
                pairCount[iBin]+=1;
            } else throw MESSAGE << "Building Cartesian Harmonic Correlations not implemented"<<ENDM_FATAL;
        }
    } else {
        for (int iBin=0; iBin<result.ndata; ++iBin ){
            double q = result.midBin( iBin );
            if (result.realpart) weight =   SpherHarmonics::ReYlm(result.l,result.m,r[0],r[1],r[2]);
            else                 weight = - SpherHarmonics::ImYlm(result.l,result.m,r[0],r[1],r[2]);
            weight *= (kernelPtr->GetValue(result.l, q, rinv));
            result.data[iBin]   += weight;
            result.uncert[iBin] += weight*weight;
            pairCount[iBin]     += 1;
        }
    }
    return true;
 }


// ------------------ postProcessPairs ------------------
bool COSCARCorrelationGenerator1d::postProcessPairs( void ){
    cout<<"Normalizing 1D source\n";
    double weight;
    result.covmtx_is_active=false;
    for (int i=0;i<result.ndata;++i) {
        if (result.data[i]!=0.0) {
            weight = 1.0/SQRTFOURPI;
            result.data[i]=result.data[i]/weight/double(pairCount[i]);
            result.uncert[i]=sqrt(abs(result.uncert[i]/double(pairCount[i])/weight/weight - result.data[i]*result.data[i])/double(pairCount[i]));
        } else result.uncert[i]=1.0;
        if (result.l==0) result.data[i]+=1.0;
    }

    if (result.l==0) {
        cout << "Cross-checking angle averaged correlation function integral...\n";
        double int_corr=0.0,int_corr2=0.0;
        for(int i=0;i<result.ndata;++i){
            int_corr  += 4.0*PI*result.midBin(i)*result.midBin(i)*result.binWidth(i)*(result.data[i]-1.0);
            int_corr2 += (result.midBin(i)*result.midBin(i)+result.binWidth(i)*result.binWidth(i)/12.0)*result.binWidth(i)*4.0*PI*(result.data[i]-1.0);
        }
        cout << "   Naive integral of correlation: " << int_corr << "\n";
        cout << "   Not so naive integral of correlation: " << int_corr2 << "\n";
        cout << "   Number of good pairs = "<<totalPairs<<"\n";
    }
    return true;
}

// ------------------ generateCorrelation ------------------
CCorrFtn1dHisto COSCARCorrelationGenerator1d::generateCorrelation( vector<COSCARLine> plist, const parameterMap& m ){

    // Initialize this instance of the correlation generator
    Read(m);
    
    // Initialize the source we are generating
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
    
    // Now build the correlation
    accumulatePairs();
    
    return result;
}
