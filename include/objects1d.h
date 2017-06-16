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
#ifndef OBJINTERFACES_1D_H
#define OBJINTERFACES_1D_H

#include "parametermap.h"

//-------------------------------------------------------------
///  Simple interface for all 1d objects
//-------------------------------------------------------------
class CObject1d {

public:
    // indices for various schemes
    int l,m;
    int lx,ly,lz;
    bool realpart;
    bool sphrHarm; //true = spherical harmonic, false = Cartesian harmonic

    // Constructors
    CObject1d(void): 
        l(0), m(0), lx(0), ly(0), lz(0), realpart(true), sphrHarm(true){}
    CObject1d(int _l, int _m, bool r): 
        l(_l), m(_m), lx(0), ly(0), lz(0), realpart(r), sphrHarm(true){}
    CObject1d(int _lx, int _ly, int _lz): 
        l(0), m(0), lx(_lx), ly(_ly), lz(_lz), realpart(true), sphrHarm(false){}
    CObject1d(const CObject1d& A): 
        l(A.l), m(A.m), lx(A.lx), ly(A.ly), lz(A.lz), realpart(A.realpart), sphrHarm(A.sphrHarm){}

    // Destructor
    virtual ~CObject1d(void){ }
    
    // Read/write to parameter map
    virtual bool Read(const parameterMap& m);
    virtual bool Write(parameterMap& m);

    // Main interface to source value (YOU MUST OVERRIDE THESE!!)
    virtual double getValue(double r) const=0;
    virtual double getError(double r) const=0;
    virtual double getCovariance(double r1, double r2) const {return 0.;} //don't make virtual yet, not implemented anywhere

    // Main interface to source value (DO NOT OVERRIDE THESE!!)
    virtual double getValue(double r){return const_cast<const CObject1d*>(this)->getValue(r);}
    virtual double getError(double r){return const_cast<const CObject1d*>(this)->getError(r);}
    virtual double getCovariance(double r1, double r2){return const_cast<const CObject1d*>(this)->getCovariance(r1,r2);}

    // CopyState
    virtual void CopyState(const CObject1d& A);

};

#endif
