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
#ifndef __OBJECTIOTEST_H
#define __OBJECTIOTEST_H

#include <string>
#include <cxxtest/TestSuite.h>
#include "parametermap.h"
#include "any.h"

using namespace std;



class ObjectIOTest : public CxxTest::TestSuite
{
public:
    void testStreamExtraction()
    {
        parameterMap com;
        ifstream f( "testData_0_0_re.coral" );
        f >> com;
        
        CCorrFtn3dSphr c;
        c.Read(com);
        
        for ( CCorrFtn3dSphr::iterator it=c.begin(); it!=c.end(); ++it) {
            cout  << it->first.l << " "  << it->first.m  << " " << boolalpha << it->first.realpart << endl;
        }     
        double r=10.;
        for (int l=0;l<=c.lmax;++l){
            for (int m=0;m<=l;++m){
                cout << l << ", " << m << ", "<<r<<": ";
                cout << c.getValueSphr(l,m,r);
                cout << " +/- ";
                cout << c.getErrorSphr(l,m,r) << endl;
            }
        }
        
        parameterMap tmp;
        c(0,0,true).Write(tmp);
        cout << tmp << endl;
 
//        TS_ASSERT_EQUALS( aMap.getItem( "aNum", 42 ), 43 );      
//        TS_ASSERT_EQUALS( aMap.getItem( "aSmallString", string( "john" ) ), "fred" );      
//        TS_ASSERT_EQUALS( parameter::getS( aMap, string( "aBigString" ) ), "insert\nlots of \ntext here!" );
    }
    
    
