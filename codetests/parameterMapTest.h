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

using namespace std;


/* THE FOLLOWING SHOULD BE PUT IN A FILE CALLED "testData_0_0_re.coral" FOR DATA TESTING PURPOSES

  string testData = "  bool bigQ  false
  bool fixed_width_bins  true
  double dx  2
  int l  0
  int m  0
  int ndata  50
  string particle1  pi+
  string particle2  pi+
  double xoffset  1
  bool realpart  true
  matrix_double covmtx_datablock {
    2.51848e-07 5.58025e-08 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    5.58025e-08 2.98916e-08 3.18493e-09 1.99752e-12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 3.18493e-09 6.44446e-09 1.62477e-09 3.87047e-12 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 1.99752e-12 1.62477e-09 1.02175e-08 3.63401e-09 1.47703e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 3.87047e-12 3.63401e-09 2.30421e-08 7.87752e-09 1.5815e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 1.47703e-11 7.87752e-09 4.12536e-08 1.24978e-08 3.58578e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 1.5815e-11 1.24978e-08 6.08753e-08 1.68683e-08 4.37008e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 3.58578e-11 1.68683e-08 7.83642e-08 2.15553e-08 5.28832e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 4.37008e-11 2.15553e-08 9.87286e-08 2.68125e-08 5.98084e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 5.28832e-11 2.68125e-08 1.2088e-07 3.25954e-08 5.91015e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 5.98084e-11 3.25954e-08 1.42573e-07 3.76012e-08 9.97478e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 5.91015e-11 3.76012e-08 1.64146e-07 4.31442e-08 9.0208e-11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 9.97478e-11 4.31442e-08 1.87215e-07 4.92182e-08 1.09469e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 9.0208e-11 4.92182e-08 2.12729e-07 5.49768e-08 1.06162e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 1.09469e-10 5.49768e-08 2.31325e-07 5.96386e-08 1.28812e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 1.06162e-10 5.96386e-08 2.49117e-07 6.37269e-08 1.50442e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.28812e-10 6.37269e-08 2.66332e-07 6.77773e-08 1.26908e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.50442e-10 6.77773e-08 2.84128e-07 7.11558e-08 1.60445e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.26908e-10 7.11558e-08 2.97197e-07 7.37682e-08 1.53659e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.60445e-10 7.37682e-08 3.03955e-07 7.53772e-08 1.6393e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.53659e-10 7.53772e-08 3.0747e-07 7.61184e-08 1.55067e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.6393e-10 7.61184e-08 3.09034e-07 7.62287e-08 1.50871e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.55067e-10 7.62287e-08 3.08131e-07 7.58923e-08 1.73925e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.50871e-10 7.58923e-08 3.05552e-07 7.44759e-08 1.55247e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.73925e-10 7.44759e-08 3.00365e-07 7.39004e-08 1.58335e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.55247e-10 7.39004e-08 2.98791e-07 7.25294e-08 1.53357e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.58335e-10 7.25294e-08 2.92711e-07 7.12442e-08 1.47619e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.53357e-10 7.12442e-08 2.86052e-07 6.99554e-08 1.55126e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.47619e-10 6.99554e-08 2.82306e-07 6.81395e-08 1.44088e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.55126e-10 6.81395e-08 2.75053e-07 6.70645e-08 1.47204e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.44088e-10 6.70645e-08 2.70828e-07 6.57455e-08 1.42386e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.47204e-10 6.57455e-08 2.6469e-07 6.46618e-08 1.34742e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.42386e-10 6.46618e-08 2.60484e-07 6.3097e-08 1.39248e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.34742e-10 6.3097e-08 2.5418e-07 6.19195e-08 1.30198e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.39248e-10 6.19195e-08 2.4995e-07 6.0785e-08 1.32749e-10 0 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.30198e-10 6.0785e-08 2.43685e-07 5.9348e-08 1.26403e-10 0 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.32749e-10 5.9348e-08 2.39529e-07 5.82837e-08 1.26903e-10 0 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.26403e-10 5.82837e-08 2.35212e-07 5.7118e-08 1.26212e-10 0 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.26903e-10 5.7118e-08 2.30317e-07 5.60967e-08 1.11708e-10 0 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.26212e-10 5.60967e-08 2.25623e-07 5.51105e-08 1.22705e-10 0 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.11708e-10 5.51105e-08 2.23343e-07 5.41184e-08 1.15702e-10 0 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.22705e-10 5.41184e-08 2.18478e-07 5.37231e-08 1.19606e-10 0 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.15702e-10 5.37231e-08 2.16995e-07 5.27625e-08 1.09706e-10 0 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.19606e-10 5.27625e-08 2.12945e-07 5.20667e-08 1.09365e-10 0 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.09706e-10 5.20667e-08 2.10625e-07 5.17345e-08 1.18333e-10 0 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.09365e-10 5.17345e-08 2.10993e-07 5.13528e-08 1.08171e-10 0 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.18333e-10 5.13528e-08 2.07395e-07 5.09635e-08 1.11061e-10 0 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.08171e-10 5.09635e-08 2.07202e-07 5.05872e-08 1.07533e-10 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.11061e-10 5.05872e-08 2.05801e-07 5.02019e-08 
    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1.07533e-10 5.02019e-08 1.82213e-07 
  }

  matrix_double ydy_datablock {
    0.612201 0.000501845 
    1.09214 0.000172892 
    1.443 8.02774e-05 
    1.57293 0.000101082 
    1.61572 0.000151796 
    1.61304 0.00020311 
    1.58132 0.000246729 
    1.53607 0.000279936 
    1.48068 0.000314211 
    1.42101 0.000347678 
    1.35824 0.000377588 
    1.2966 0.000405149 
    1.24022 0.000432684 
    1.19088 0.000461226 
    1.14732 0.000480963 
    1.10959 0.000499116 
    1.08127 0.000516074 
    1.05662 0.000533036 
    1.03882 0.000545158 
    1.02671 0.000551321 
    1.01715 0.0005545 
    1.0101 0.000555909 
    1.00659 0.000555096 
    1.00309 0.000552768 
    1.00136 0.000548055 
    1.00052 0.000546618 
    0.997711 0.000541028 
    0.996951 0.000534838 
    0.997076 0.000531325 
    0.997047 0.000524455 
    0.997869 0.000520411 
    0.997582 0.00051448 
    0.997733 0.000510376 
    0.998225 0.000504162 
    0.99841 0.00049995 
    0.998412 0.000493645 
    0.99759 0.000489417 
    0.997933 0.000484987 
    0.99718 0.000479914 
    0.996496 0.000474998 
    0.997739 0.000472591 
    0.997875 0.000467417 
    0.998687 0.000465827 
    0.997523 0.00046146 
    0.998127 0.00045894 
    0.999743 0.000459339 
    0.998683 0.000455406 
    0.999538 0.000455194 
    0.998244 0.000453653 
    0.999475 0.000426864 
  }
"
*/


class ParameterMapTest : public CxxTest::TestSuite
{
public:
    void testSetAndGetItemDefault()
    {
        // This is the easiest way to use parameter maps, use the set function in the namespace (set can get the type 
        // of the value from the value itself since it gets stuck in an any object).  The getItem() member function
        // knows what type the any points to since it has a default of the same type.
        parameterMap aMap;
        parameter::set( aMap, "aNum", 43 ); // note, this stores 43 as an int!
        parameter::set( aMap, "aSmallString", "fred" );
        parameter::set( aMap, "aBigString", "insert\nlots of \ntext here!" );        
        TS_ASSERT_EQUALS( aMap.getItem( "aNum", 42 ), 43 );      
        TS_ASSERT_EQUALS( aMap.getItem( "aSmallString", string( "john" ) ), "fred" );      
        TS_ASSERT_EQUALS( parameter::getS( aMap, string( "aBigString" ) ), "insert\nlots of \ntext here!" );
    }
    
    void testSetAndGetItem()
    {
        // This is the 2nd easiest way to use parameter maps, use the set function in the namespace.  
        // The getItem() member function takes a template typename so we know which getItem() to use.
        parameterMap aMap;
        parameter::set( aMap, "aNum", 43.0 ); // note, here 43 is a double@
        parameter::set( aMap, "aSmallString", string( "fred" ) );
        parameter::set( aMap, "aBigString", string( "insert\nlots of \ntext here!" ) );
        TS_ASSERT_EQUALS( aMap.getItem< double >( string( "aNum" ) ), 43.0 );      
        TS_ASSERT_EQUALS( aMap.getItem< string >( string( "aBigString" ) ), "insert\nlots of \ntext here!" );      
    }

    void testSTLSetAndGet()
    {
        // This is yet another way to use the parameter map.  Instead of calling a set function, you can assign
        // directly using the operator[] to get a reference since the parametermap inherits from the std::map class.  
        // You can pull the data out by using operator[] function as well, but you get back an any, which you have 
        // to cast into the type you really need.  Also, any_cast'ing is very picky -- you can't any_cast an int into a double.
        parameterMap aMap;
        aMap[ "aNum" ] = 43 ;
        aMap[ "aSmallString" ] = string( "fred" );
        aMap[ "aBigString" ] = string( "insert\nlots of \ntext here!" );
        TS_ASSERT_EQUALS( boost::any_cast< int >( aMap[ "aNum" ] ), 43 );      
        TS_ASSERT_EQUALS( boost::any_cast< string >( aMap[ "aBigString" ] ), "insert\nlots of \ntext here!" );      
    }

    void testSetAndGet()
    {
        // Here is another way to do the gets: there is one getX function for each type (getD gets a double, etc).
        // getX() knows to call getItem< X >().
        parameterMap aMap;
        parameter::set( aMap, "aNum", 43.0 );
        parameter::set( aMap, "aSmallString", "fred" );
        parameter::set( aMap, "aBigString", "insert\nlots of \ntext here!" );        
        TS_ASSERT_EQUALS( parameter::getD( aMap, "aNum"), 43.0 );
        TS_ASSERT_EQUALS( parameter::getS( aMap, "aSmallString" ), "fred" );
        TS_ASSERT_EQUALS( parameter::getS( aMap, "aBigString" ), "insert\nlots of \ntext here!" );
    }
    
    void testSetAndGetWithDefaults()
    {
        // Same as above, but with default arguments.
        parameterMap aMap;
        parameter::set( aMap, "aNum", 43 );
        parameter::set( aMap, "aSmallString", "fred" );
        parameter::set( aMap, "aBigString", "insert\nlots of \ntext here!" );        
        TS_ASSERT_EQUALS( parameter::getI( aMap, "aNum", 43 ), 43 );
        TS_ASSERT_EQUALS( parameter::getS( aMap, "anotherSmallString", "fred" ), "fred" );
        TS_ASSERT_EQUALS( parameter::getS( aMap, "aBigString", "fred" ), "insert\nlots of \ntext here!" );
    }

    
    void testFileReading()
    {
        // This time we read from a map using the stream extraction operator (operator>>).  Just to check 
        // that the stream extraction is behaving properly, we getX a few elements of the map.  Note, 
        // we use the getM function to get a matrix of doubles.
        string datafile( "testData_0_0_re.coral" );
        ifstream f( datafile.c_str() );
        parameterMap aMap; 
        if ( !f.good() ) cerr << "\nOpening data file " << datafile << " failed!" << endl;
        while ( f.good() ) f >> aMap;
        f.close();
        TS_ASSERT_EQUALS( parameter::getB( aMap, "fixed_width_bins" ), true );
        TS_ASSERT_EQUALS( parameter::getS( aMap, "particle2" ), "pi+" );
        TS_ASSERT_EQUALS( parameter::getI( aMap, "l" ), 0 );
        TS_ASSERT_EQUALS( parameter::getM( aMap, "ydy_datablock" )[3][1], 0.000101082 );
    }

};


#endif 
