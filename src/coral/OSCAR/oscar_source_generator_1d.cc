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
#include "oscar_source_generator_1d.h"
#include "constants.h"
#include "tnt_array1d_utils.h"
#include "misc.h"
#include "message.h"
#include "sf.h"


// ------------------ addOnePair ------------------
bool COSCARSourceGenerator1d::addOnePair( const COSCARLine& p1, const COSCARLine& p2 ){
    // This is a good pair, let's bin it up the source function(s)...
    int iBin = result.whatBin(rinv,false);
    if (iBin>=0 && iBin<result.ndata) {
        if ( result.sphrHarm ) {
            double r[3];
            double weight;
            getSideOutLong( r_cm, r );
            double rinv2=rinv*rinv; //r_cm[0]*r_cm[0]+r_cm[1]*r_cm[1]+r_cm[2]*r_cm[2];
            double binWidth = result.binWidth( iBin );
            if (rinv < 1e-6) {
                rinv2 = result.midBin( iBin )*result.midBin( iBin );
                MESSAGE << "Bad pair, set rinv^2 to " << rinv2 << ENDM_WARN; 
            }
            if (result.realpart) weight =   SpherHarmonics::ReYlm(result.l,result.m,r[0],r[1],r[2])/rinv2/binWidth;
            else                 weight = - SpherHarmonics::ImYlm(result.l,result.m,r[0],r[1],r[2])/rinv2/binWidth;
            result.data[iBin]   += weight;
            result.uncert[iBin] += weight*weight;
        } else throw MESSAGE << "Building Cartesian Harmonic Sources not implemented"<<ENDM_SEVERE;
    }
    return true;
 }


// ------------------ postProcessPairs ------------------
bool COSCARSourceGenerator1d::postProcessPairs( void ){
    cout<<"    Normalizing 1D source\n";
    result.covmtx_is_active=false;
    for (int i=0;i<result.ndata;++i) result.data[i]=result.data[i]/SQRTFOURPI/double(totalPairs);
    for (int i=0;i<result.ndata;++i) {
        double u = result.uncert[i]/double(totalPairs)/4.0/PI-result.data[i]*result.data[i];
        if ( u < -1e-6 ) throw MESSAGE <<"Bad u: "<<u<<ENDM_FATAL;
        result.uncert[i]=sqrt(abs(u/double(totalPairs)));
    }
    double dintegral=0.0, dintegral2=0.0;
    if (result.l==0) {
        cout << "    Cross-checking angle averaged source function integral...\n";
        double int_sou=0.0,int_sou2=0.0;
        for(int i=0;i<result.ndata;++i){
            dintegral += pow(result.uncert[i]*4.0*PI*result.midBin(i)*result.midBin(i)*result.binWidth(i),2.0);
            dintegral2 += pow(result.uncert[i]*4.0*PI*(result.midBin(i)*result.midBin(i)+result.binWidth(i)*result.binWidth(i)/12.0),2.0);
            int_sou  += 4.0*PI*result.midBin(i)*result.midBin(i)*result.binWidth(i)*result.data[i];
            int_sou2 += (result.midBin(i)*result.midBin(i)+result.binWidth(i)*result.binWidth(i)/12.0)*result.binWidth(i)*4.0*PI*result.data[i];
        }
        cout << "        Naive integral of source: " << int_sou << " +/- " << dintegral << "\n";
        cout << "        Not so naive integral of source: " << int_sou2 << " +/- " << dintegral2 << "\n";
        cout << "        Number of good pairs = "<<totalPairs<<"\n";
    }
    return true;
}

// ------------------ generateSource ------------------
CSourceFtn1dHisto COSCARSourceGenerator1d::generateSource( vector<COSCARLine> plist, const parameterMap& m ){

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
