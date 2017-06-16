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
#ifndef __OSCAR_CORRELATION_GEN3DCART_H__
#define __OSCAR_CORRELATION_GEN3DCART_H__

#include "parametermap.h"
#include "oscar_correlation_generatorbase.h"
#include "corr3d_histo.h"

/** 
  * This class handles the computation of a three-dimensional correlation function in the pair center-of-mass frame:
  * \f$C(\vec{q}) = 1+\int d^3r' K(\vec{q},\vec{r'}) \int dt' \int d^4R D(R+r/2, \vec{P}) D(R-r/2, \vec{P})\f$
  *
  */
class COSCARCorrelationGenerator3dCart: public COSCARCorrelationGeneratorBase{

    public:
        
        COSCARCorrelationGenerator3dCart( void ): COSCARCorrelationGeneratorBase(), result(), 
            maxS(0.0), minS(0.0), maxO(0.0), minO(0.0), maxL(0.0), minL(0.0){for (int i=0;i<3;++i) fold[i]=false;}
        ~COSCARCorrelationGenerator3dCart( void ){}

        bool Read( const parameterMap& m );
        bool Write( parameterMap& m);
        
        bool addOnePair( const COSCARLine& p1, const COSCARLine& p2 );
        bool postProcessPairs( void );

        CCorrFtn3dHisto generateCorrelation( vector<COSCARLine> plist, const parameterMap& m );       

        // Define tmp variables 
        CCorrFtn3dHisto result;
        double maxS,minS,maxO,minO,maxL,minL;
        bool fold[3];
        TNT::Array1D<int> pairCount;
};

#endif
