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
#ifndef __OSCAR_CORRELATION_GEN1D_H__
#define __OSCAR_CORRELATION_GEN1D_H__

#include "parametermap.h"
#include "oscar_correlation_generatorbase.h"
#include "corr1d_histo.h"
#include "tnt_array1d.h"

/** 
  * This class handles the computation of a one-dimensional correlation function, in the pair center-of-mass frame, 
  * projected onto spherical or Cartesian harmonics:
  * \f$C_{\ell m}(|\vec{q}|) = \delta_{\ell 0} + \sqrt{4\pi} \int d^3 r' K_\ell (|\vec{q}|, |\vec{r}'|) Y^*_{\ell m}(\Omega_{\vec{r}'})\int dt' \int d^4R D(R+r/2, \vec{P}) D(R-r/2, \vec{P})\f$
  *
  */
class COSCARCorrelationGenerator1d: public COSCARCorrelationGeneratorBase{

    public:
        
        COSCARCorrelationGenerator1d( void ): COSCARCorrelationGeneratorBase(), result(), pairCount(0){}
        ~COSCARCorrelationGenerator1d( void ){}

        bool addOnePair( const COSCARLine& p1, const COSCARLine& p2 );
        bool postProcessPairs( void );

        CCorrFtn1dHisto generateCorrelation( vector<COSCARLine> plist, const parameterMap& m );       
        
        // Define tmp variables 
        CCorrFtn1dHisto result;
        TNT::Array1D<int> pairCount;
};

#endif
