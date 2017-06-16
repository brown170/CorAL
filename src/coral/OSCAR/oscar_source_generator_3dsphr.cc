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
#include "oscar_source_generator_3dsphr.h"
#include "constants.h"
#include "tnt_array1d_utils.h"
#include "misc.h"
#include "message.h"
#include "sf.h"

// ------------------ Read ------------------
bool COSCARSourceGenerator3dSphr::Read( const parameterMap& m ){
    lmax = parameter::getI(m,"lmax",lmax);
    return COSCARSourceGeneratorBase::Read(m);
}


// ------------------ Write ------------------
bool COSCARSourceGenerator3dSphr::Write( parameterMap& m){
    parameter::set(m,"lmax",lmax);
    return COSCARSourceGeneratorBase::Write(m);
}


// ------------------ addOnePair ------------------
bool COSCARSourceGenerator3dSphr::addOnePair( const COSCARLine& p1, const COSCARLine& p2 ){
    // bin up 1D source
    double r[3];
    int iBin; 
    getSideOutLong( r_cm, r );
    double rinv2=rinv*rinv;//r_cm[0]*r_cm[0]+r_cm[1]*r_cm[1]+r_cm[2]*r_cm[2];
    for ( CSourceFtn3dSphr<CSourceFtn1dHisto>::iterator it = result.begin(); it!=result.end(); ++it ){
        iBin = it->second.whatBin(rinv,false);
        if ( rinv < 1e-6 ) {
            rinv2 = it->second.midBin( iBin );
            rinv2 *= rinv2;
            MESSAGE << "Bad pair, set rinv^2 to " << rinv2 << ENDM_WARN; 
        }
        if (iBin>=0 && iBin<it->second.ndata) {
            if ( it->second.sphrHarm ) {
                double weight;
                if (it->second.realpart) weight =   SpherHarmonics::ReYlm(it->second.l,it->second.m,r[0],r[1],r[2])/rinv2;
                else                     weight = - SpherHarmonics::ImYlm(it->second.l,it->second.m,r[0],r[1],r[2])/rinv2;
                it->second.data[iBin]   += weight;
                it->second.uncert[iBin] += (weight*weight);
            } else throw MESSAGE << "Building Cartesian Harmonic Sources for l>0 not implemented"<<ENDM_FATAL;
        }
    }
    return true;
}


// ------------------ postProcessPairs ------------------
bool COSCARSourceGenerator3dSphr::postProcessPairs( void ){

    cout<<"    Normalizing 1D sources\n";

    for ( CSourceFtn3dSphr<CSourceFtn1dHisto>::iterator it = result.begin(); it!=result.end(); ++it ){
        it->second.covmtx_is_active=false;
        double rweight,rval,dr;
        for (int i=0;i<it->second.ndata;++i) {
            rval    = it->second.midBin(i);
            dr      = it->second.binWidth(i);
            //rweight = (rval*rval+dr*dr/12.0)*dr*SQRTFOURPI;
            rweight = dr*SQRTFOURPI;
            it->second.data[i]   = (it->second.data[i])/rweight/double(totalPairs);
            double u = it->second.uncert[i] = it->second.uncert[i]/double(totalPairs)/rweight/rweight-it->second.data[i]*it->second.data[i];
            if ( u < -1e-6 ) throw MESSAGE <<"Bad u: "<<u<<ENDM_FATAL;
            it->second.uncert[i] = sqrt(abs(u/double(totalPairs)));
        }

        if (it->second.l==0) {
            double dintegral=0.0, dintegral2=0.0;
            cout << "    Cross-checking angle averaged source function integral...\n";
            double int_sou=0.0,int_sou2=0.0;
            for(int i=0;i<it->second.ndata;++i){
                dintegral += pow(it->second.uncert[i]*4.0*PI*it->second.midBin(i)*it->second.midBin(i)*it->second.binWidth(i),2.0);
                dintegral2 += pow(it->second.uncert[i]*4.0*PI*(it->second.midBin(i)*it->second.midBin(i)+it->second.binWidth(i)*it->second.binWidth(i)/12.0),2.0);
                int_sou  += 4.0*PI*(it->second.midBin(i))*(it->second.midBin(i))*(it->second.binWidth(i))*(it->second.data[i]);
                int_sou2 += (it->second.midBin(i)*(it->second.midBin(i))+it->second.binWidth(i)*(it->second.binWidth(i))/12.0)*(it->second.binWidth(i))*4.0*PI*(it->second.data[i]);
            }
            cout << "        Naive integral of source: " << int_sou << " +/- " << dintegral << "\n";
            cout << "        Not so naive integral of source: " << int_sou2 << " +/- " << dintegral2 << "\n";
            cout << "        Number of good pairs = "<<totalPairs<<"\n";
        }
    }
    return true;
}

// ------------------ generateSource ------------------
CSourceFtn3dSphr<CSourceFtn1dHisto> COSCARSourceGenerator3dSphr::generateSource( vector<COSCARLine> plist, const parameterMap& m ){

    // Initialize the source we are generating
    parameterMap souMap; 
    souMap = parameter::getMap( m,"source_settings" );
    result.Read(souMap);
    result.fillTerms();
    for ( CSourceFtn3dSphr< CSourceFtn1dHisto >::iterator it=result.begin(); it!=result.end(); ++it ) it->second.Read(souMap);
    
    // Initialize this instance of the source generator
    Read(m);

    // Filter the pairs
    filterParticles( plist );

    // Now build the source
    accumulatePairs();

    return result;
}
