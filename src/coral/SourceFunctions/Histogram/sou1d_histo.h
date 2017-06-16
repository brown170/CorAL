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
#ifndef SOURCES_1DHISTO_H
#define SOURCES_1DHISTO_H

#include "soubase.h"
#include "histogram1d.h"

class CSourceFtn1dHisto: public CSourceFtnBase, public CHistogram1d {
public:
    // create/destroy
    CSourceFtn1dHisto(string p1="", string p2="", 
        int l=0, int m=0, bool repart=true,int n=0): 
        CSourceFtnBase(p1,p2), CHistogram1d(l,m,repart,n) {}
    CSourceFtn1dHisto(const CSourceFtn1dHisto& A): CSourceFtnBase(A), CHistogram1d(A){}
            
    // read/write to parameter map
    bool Read(const parameterMap& s){return CSourceFtnBase::Read(s)&&CHistogram1d::Read(s);}
    bool Write(parameterMap& s){return CSourceFtnBase::Write(s)&&CHistogram1d::Write(s);}
        
    // copy
    void CopyState(const CSourceFtn1dHisto& A) 
        {CSourceFtnBase::CopyState(A); CHistogram1d::CopyState(A);}
            
};

#endif
