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
#ifndef __PARAMETERMAPTEST_H
#define __PARAMETERMAPTEST_H

#include <string>
#include <cxxtest/TestSuite.h>
#include "parametermap.h"
#include "any.h"
#include "corr1d_histo.h"
#include "imager.h"

using namespace std;


class ParameterMapTest : public CxxTest::TestSuite
{
public:
    void testReadCorrelation(void){
        if ( !file_exists( "testData_0_0_re.coral") ) throw MESSAGE << "File "<< "testData_0_0_re.coral" <<" not found. Halting!!"<<ENDM_FATAL;
        CCorrFtn1dHisto corrin( "pi0", "pi0", false );
        parameterMap corrin_params;
        ifstream corrin_file( "testData_0_0_re.coral" );
        corrin_file >> corrin_params;
        corrin.Read( corrin_params );
        TS_ASSERT_EQUALS( parameter::getB( corrin_params, "fixed_width_bins" ), true );
        TS_ASSERT_EQUALS( parameter::getS( corrin_params, "particle2" ), "pi+" );
        TS_ASSERT_EQUALS( parameter::getI( corrin_params, "l" ), 0 );
        TS_ASSERT_EQUALS( parameter::getM( corrin_params, "ydy_datablock" )[3][1], 0.000101082 );
        TS_ASSERT_EQUALS( corrin.covmtx[3][3], 1.02175e-08 );
        TS_ASSERT_EQUALS( corrin.particle1, "pi+" );
        TS_ASSERT_EQUALS( corrin.particle2, "pi+" );
        TS_ASSERT_EQUALS( corrin.bigQ, false );
        TS_ASSERT_EQUALS( corrin.fixed_width_bins, true );
        TS_ASSERT_EQUALS( corrin.ndata, 50 );       
    }
    
    void testImage(void){
        // Read in and initialize the correlation function
        if ( !file_exists( "testData_0_0_re.coral") ) throw MESSAGE << "File "<< "testData_0_0_re.coral" <<" not found. Halting!!"<<ENDM_FATAL;
        CCorrFtn1dHisto corrin( "pi0", "pi0", false );
        parameterMap corrin_params;
        ifstream corrin_file( "testData_0_0_re.coral" );
        corrin_file >> corrin_params;
        corrin.Read( corrin_params );
        
        // Set the parameters for the imager
        parameterMap image_com;
        image_com[ "constrain_origin" ] = true;
        image_com[ "constrain_rmax_zero" ] = false;
        image_com[ "qmax" ] = 70.0;
        image_com[ "rmax" ] = 60.0; 
        image_com[ "numcoeffs" ] = 10;
        image_com[ "hbt_only" ] = false;

        // Now actually image
        CSourceFtn1dLegendrePoly souout( image_le( corrin, image_com ) );
        
        // Save imaged source
        parameterMap sou_map;
        souout.Write( sou_map );
        ofstream sourceOut( "pion_legendre_source.dat" );
        sourceOut << sou_map;
        
        // Something plottable
        ofstream dataOut( "legendre_data.dat" );
        for (int i=0; i<52; ++i){
            double x = (double)i;
            dataOut << x << "  " << souout.getValue(x) << "  " << souout.getError(x)<<endl;
        }
    }
};


#endif 
