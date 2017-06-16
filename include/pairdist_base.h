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
#ifndef PAIRDIST_BASE_H
#define PAIRDIST_BASE_H

#include <string>
#include "parametermap.h"

using namespace std;

//------------------------------------------------------------------------
//  Base class for all pair momentum distributions (yes it is a copy of
//  correlation classes)
//------------------------------------------------------------------------
class CPairDistributionBase {

public:
    // public member data
    string particle1;
    string particle2;
    bool bigQ;
    
    // create
    CPairDistributionBase(string p1 = "", string p2 = "", bool bQ = false) :
        particle1(p1), particle2(p2), bigQ(bQ) {}
    CPairDistributionBase(const CPairDistributionBase& m) : 
        particle1(m.particle1), particle2(m.particle2), bigQ(m.bigQ){}
    virtual ~CPairDistributionBase(void){}


    // read/write to parameterMap map
    virtual bool Read(const parameterMap& s);
    virtual bool Write(parameterMap& s);

    bool likepair(void) const{return particle1 == particle2;}

    // CopyState
    virtual void CopyState(const CPairDistributionBase& A)
        {particle1=A.particle1; particle2=A.particle2; bigQ=A.bigQ;}

};

#endif
