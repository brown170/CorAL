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
#ifndef OBJINTERFACES_3D_H
#define OBJINTERFACES_3D_H

#include "parametermap.h"

//-------------------------------------------------------------
///  Simple interface for all 3d objects
//-------------------------------------------------------------
class CObject3d  {

public:
    // Constructors
    CObject3d(void){}
    CObject3d(const CObject3d& A){}

    // Destructor
    virtual ~CObject3d(void){ }
    
    // Read/write to parameter map
    virtual bool Read(const parameterMap& m){return true;}
    virtual bool Write(parameterMap& m){return true;}

    // Main interface to source value (YOU SHOULD OVERRIDE THESE!!)
    virtual double getValueCart(double x, double y, double z) const = 0;
    virtual double getValueSphr(double r, double theta, double phi) const = 0;
    virtual double getErrorCart(double x, double y, double z) const = 0;
    virtual double getErrorSphr(double r, double theta, double phi) const = 0;

    // Main interface to source value (DO NOT OVERRIDE THESE!!)
    virtual double getValueCart(double x, double y, double z) {return const_cast<const CObject3d*>(this)->getValueCart(x, y, z);}
    virtual double getValueSphr(double r, double theta, double phi) {return const_cast<const CObject3d*>(this)->getValueCart(r, theta, phi);}
    virtual double getErrorCart(double x, double y, double z) {return const_cast<const CObject3d*>(this)->getErrorCart(x, y, z);}
    virtual double getErrorSphr(double r, double theta, double phi) {return const_cast<const CObject3d*>(this)->getErrorCart(r, theta, phi);}

    // CopyState
    virtual void CopyState(const CObject3d& A){}

    // readTerms 
    //! Read terms from disk, must override as 
    //! naming scheme dependents on term type 
    //! (provided derived class even has terms!)
    virtual bool readTerms( void ){ return true; }

    // writeTerms 
    //! Write terms to disk, must override as 
    //! naming scheme dependents on term type 
    //! (provided derived class even has terms!)
    virtual bool writeTerms( void ){ return true; }

};

#endif
