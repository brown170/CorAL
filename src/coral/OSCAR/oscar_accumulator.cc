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
#include "oscar_accumulator.h"
#include "misc.h"

// ------------------ Read ------------------
bool COSCARAccumulator::Read( const parameterMap& m ){
    pid1 = parameter::getI( m, "pid1", 211 );
    pid2 = parameter::getI( m, "pid2", 211 );
    max_number_pairs = parameter::getI( m, "max_number_pairs", 1000000 );
    // Determine which of Side, Out and Long is x, y, and z
    iSide = parameter::getI(m, "cartesian_side_index", iSide);
    iOut  = parameter::getI(m, "cartesian_out_index",  iOut);
    iLong = parameter::getI(m, "cartesian_long_index", iLong);
    return true;
}

// ------------------ Write ------------------
bool COSCARAccumulator::Write( parameterMap& m){
    parameter::set( m, "pid1", pid1 );
    parameter::set( m, "pid2", pid2 );
    parameter::set( m, "max_number_pairs", max_number_pairs );
    parameter::set(m,"cartesian_side_index",iSide);
    parameter::set(m,"cartesian_out_index", iOut);
    parameter::set(m,"cartesian_long_index",iLong);
    return true;
}

// ------------------ filterParticles ------------------
void COSCARAccumulator::filterParticles( vector<COSCARLine> plist ){
    cout << "    Filtering out all particles that do not have PID = "<<pid1;
    if (pid1 == pid2) {
        particleList1 = filterPID( plist, pid1 );
        particleList2.clear();
    } else {
        cout << " or PID = "<<pid2;
        particleList1 = filterPID( plist, pid1 );
        particleList2 = filterPID( plist, pid2 );
    }
    cout << "     ... done."<<endl;
    cout << "        There are "<<particleList1.size()<<" particles in list 1 and "
          <<particleList2.size()<<" particles in list 2. "<<endl;
    return;
}

// ------------------ accumulatePairs ------------------
void COSCARAccumulator::accumulatePairs( void ){
    totalPairs = 0;
    // loop over pairs
    for (vector<COSCARLine>::iterator p1=particleList1.begin(); p1!=particleList1.end(); ++p1){
        vector<COSCARLine>::iterator p2, p2_begin, p2_end;
        if (pid1==pid2) {
            p2_begin = p1;++p2_begin;//particleList1.begin(); //avoids double counting
            p2_end   = particleList1.end();
        } else {
            p2_begin = particleList2.begin();
            p2_end   = particleList2.end();
        }
        p2 = p2_begin;
        while ( (totalPairs<max_number_pairs) && (p2!=p2_end) ) {
            // Check if pair is good 
            this->setCOMVariables( *p1, *p2 );
            if ( this->pairIsGood( *p1, *p2 ) ){  
                totalPairs+=1;
                this->addOnePair( *p1, *p2 );
            }
            ++p2;
        }
        if (totalPairs>=max_number_pairs) break;
    }
    cout << "    totalPairs: " << totalPairs << " out of " << max_number_pairs<<" allowed"<<endl;
    this->postProcessPairs();
}

// ------------------ likePair ------------------
bool COSCARAccumulator::likePair( void ){ return pid1==pid2; }

// ------------------ getSideOutLong ------------------
void COSCARAccumulator::getSideOutLong( double* inVec, double* outVec ){ 
    // define unit three-vectors of Bertsch-Pratt coords.
    // O is unit vector || to the part of P that is perp to L 
    O[0]=P_lab[1];O[1]=P_lab[2];O[2]=0.0; // P is a 4-vector, O is a 3-vector
    double Omag = sqrt(O[0]*O[0]+O[1]*O[1]);
    O[0]/=Omag;O[1]/=Omag;
    // S is perp to both O and L
    S[0]=O[1];S[1]=-O[0];S[2]=0.0;           
    // break sep_cm into SOL coords and figure out 
    // the bins they correspond too
    double rS,rO,rL;
    rS=0.0;rO=0.0;rL=0.0;
    for (int i=0;i<3;++i) {
        rS+=inVec[i+1]*S[i];
        rO+=inVec[i+1]*O[i];
        rL+=inVec[i+1]*L[i];
    }
    // Set x, y, z according to specified coordinate system
    outVec[iSide] = rS;
    outVec[iOut]  = rO;
    outVec[iLong] = rL;
}

// ------------------ setCOMVariables ------------------
void COSCARAccumulator::setCOMVariables( const COSCARLine& p1, const COSCARLine& p2 ){
    double PSquared=0.0;
    double g[4]={1.,-1.,-1.,-1.};
    for (int i=0;i<4;++i){
        q_lab[i]   = 0.5*(p1.p[i]-p2.p[i]);  // relative mom. is 1/2 the mom. diff. of pair
        P_lab[i]   = p1.p[i]+p2.p[i];        // total mom. is total of pair mom.
        r_lab[i]   = (p1.x[i]-p2.x[i]);
        r_cm[i]    = 0.0;
        q_cm[i]    = 0.0;
        PSquared   += P_lab[i]*P_lab[i]*g[i];
    }
    for(int i=0;i<4;++i){beta[i] = g[i]*P_lab[i]/sqrt(PSquared);}
    Misc::lorentz( beta, q_lab, q_cm );
    Misc::lorentz( beta, r_lab, r_cm );
    // Loop over space-part of variable
    sqr_q_cm = q_cm[0]*q_cm[0]*g[0];
    rinv     = 0.0;
    for (int i=1; i<4; ++i) {
        sqr_q_cm += q_cm[i]*q_cm[i]*g[i];
        rinv     += r_cm[i]*r_cm[i];
    }
    rinv = sqrt( abs(rinv) );
    qinv = sqrt( abs(sqr_q_cm) );
}
