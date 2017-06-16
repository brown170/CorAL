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
#ifndef __OSCAR_SOURCE_GEN3DSPHR_H__
#define __OSCAR_SOURCE_GEN3DSPHR_H__

#include "parametermap.h"
#include "oscar_source_generatorbase.h"
#include "sou1d_histo.h"
#include "sou3d_ylm.h"

/** 
  * This class handles the computation of a three-dimensional source function, in the pair center-of-mass frame:
  * \f$S(\vec{r}') = \int dt' \int d^4R D(R+r/2, \vec{P}) D(R-r/2, \vec{P})\f$.
  * However, this class does it by first expanding the source in spherical harmonics, 
  * \f$S(\vec{r}') = \sqrt{4\pi} \sum_{\ell m} Y_{\ell m}(\Omega_{\vec{r}'}) S_{\ell m}(|\vec{r}'|)\f$
  * then solving for each term in the expansion:
  * \f$S_{\ell m}(|\vec{r}'|) = \frac{1}{\sqrt{4\pi}} \int d\Omega_{\vec{r}'} Y^*_{\ell m}(\Omega_{\vec{r}'})\int dt' \int d^4R D(R+r/2, \vec{P}) D(R-r/2, \vec{P})\f$
  *
  */
class COSCARSourceGenerator3dSphr: public COSCARSourceGeneratorBase{

    public:
        
        COSCARSourceGenerator3dSphr( void ): COSCARSourceGeneratorBase(), result(), lmax(0){}
        ~COSCARSourceGenerator3dSphr( void ){}

        bool Read( const parameterMap& m );
        bool Write( parameterMap& m);
        
        bool addOnePair( const COSCARLine& p1, const COSCARLine& p2 );
        bool postProcessPairs( void );

        CSourceFtn3dSphr<CSourceFtn1dHisto> generateSource( vector<COSCARLine> plist, const parameterMap& m );       
        
    private:
        // Define tmp variables 
        CSourceFtn3dSphr<CSourceFtn1dHisto> result;
        int lmax;
};

#endif
