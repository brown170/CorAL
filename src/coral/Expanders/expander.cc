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
#include "expander.h"
#include "cheezyparser.h"

//--------------------------------------------
// Main converter functions!  
//--------------------------------------------

CCorrFtn3dSphr expand( const CCorrFtn3dHisto& corrin, const parameterMap m ){
    // create the target object & initialize the common baseclass data
    CCorrFtn3dSphr corrout;
    corrout.CObject3d::CopyState(corrin);
    corrout.CCorrFtnBase::CopyState(corrin);
    // expand it
    CCart2SphrExpander<CCorrFtn1dHisto> ex;
    ex.Read(m);
    if (!ex.convert(corrin,corrout,corrin.likepair())) exit(-1);
    // initialize base-class data of ylm terms
    for (CCorrFtn3dSphr::iterator it=corrout.begin(); it!=corrout.end(); ++it){
        it->second.CCorrFtnBase::CopyState(corrin);
    }
    return corrout;
}

CCorrFtn3dHisto recombine( const CCorrFtn3dSphr& corrin, const parameterMap m ){
    // create the target object & initialize the common baseclass data
    CCorrFtn3dHisto corrout;
    corrout.CObject3d::CopyState(corrin);
    corrout.CCorrFtnBase::CopyState(corrin);
    // expand it
    CCart2SphrExpander<CCorrFtn1dHisto> ex;
    ex.Read(m);
    if (!ex.convert(corrin,corrout)) exit(-1);
    return corrout;
}

CSourceFtn3dSphr<CSourceFtn1dHisto> expand( const CSourceFtn3dHisto& souin, const parameterMap m ){
    // create the target object & initialize the common baseclass data
    CSourceFtn3dSphr<CSourceFtn1dHisto> souout;
    souout.CObject3d::CopyState(souin);
    souout.CSourceFtnBase::CopyState(souin);
    // expand it
    CCart2SphrExpander<CSourceFtn1dHisto> ex;
    ex.Read(m);
    if (!ex.convert(souin,souout,souin.likepair())) exit(-1);
    // initialize base-class data of ylm terms
    for (CSourceFtn3dSphr<CSourceFtn1dHisto>::iterator it=souout.begin(); it!=souout.end(); ++it){
        it->second.CSourceFtnBase::CopyState(souin);
    }
    return souout;
}

CSourceFtn3dHisto recombine( const CSourceFtn3dSphr<CSourceFtn1dHisto>& souin, const parameterMap m ){
    // create the target object & initialize the common baseclass data
    CSourceFtn3dHisto souout;
    souout.CObject3d::CopyState(souin);
    souout.CSourceFtnBase::CopyState(souin);
    // expand it
    CCart2SphrExpander<CSourceFtn1dHisto> ex;
    ex.Read(m);
    if (!ex.convert(souin,souout)) exit(-1); 
    return souout;
}

