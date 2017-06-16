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
#ifndef CORR1DHISTO_H
#define CORR1DHISTO_H

#include "corrbase.h"
#include "histogram1d.h"

class CCorrFtn1dHisto: public CCorrFtnBase, public CHistogram1d {
public:
    // create/destroy
    CCorrFtn1dHisto(string p1="", string p2="", bool bQ=false, 
        int l=0, int m=0, bool repart=true,int n=0): 
        CCorrFtnBase(p1,p2,bQ), CHistogram1d(l,m,repart,n) {}
    CCorrFtn1dHisto(const CCorrFtn1dHisto& A): CCorrFtnBase(A), CHistogram1d(A){}
            
    // read/write to parameter map
    bool Read(const parameterMap& s){return CCorrFtnBase::Read(s)&&CHistogram1d::Read(s);}
    bool Write(parameterMap& s){return CCorrFtnBase::Write(s)&&CHistogram1d::Write(s);}
        
    // copy
    void CopyState(const CCorrFtn1dHisto& A) 
        {CCorrFtnBase::CopyState(A); CHistogram1d::CopyState(A);}
            
};

#endif
