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
#ifndef __OSCAR_ACCUMULATOR_H__
#define __OSCAR_ACCUMULATOR_H__

#include "parametermap.h"
#include "oscar.h"

class COSCARAccumulator{

    public:
        
        COSCARAccumulator( void ): pid1(211), pid2(211), 
            max_number_pairs(0), totalPairs(0), particleList1(0), particleList2(0), 
            iSide(0), iOut(1), iLong(2), sqr_q_cm(0.), rinv(0.), qinv(0.0)
            {L[0]=0.0;L[1]=0.0;L[2]=1.0;}
        virtual ~COSCARAccumulator( void ){}

        virtual bool Read( const parameterMap& m );
        virtual bool Write( parameterMap& m);
        
        void filterParticles( vector<COSCARLine> plist );
        void accumulatePairs( void );
        virtual bool pairIsGood( const COSCARLine& p1, const COSCARLine& p2 )=0;
        virtual bool addOnePair( const COSCARLine& p1, const COSCARLine& p2 )=0;
        virtual bool postProcessPairs( void )=0;
        bool likePair( void );
        
        void setCOMVariables( const COSCARLine& p1, const COSCARLine& p2 );
        void getSideOutLong( double* inVec, double* outVec ); 
        
        int pid1, pid2, max_number_pairs, totalPairs;
        vector< COSCARLine > particleList1, particleList2;

    protected:
        // Define tmp variables 
        int iSide;          //! In a 3-vector, this is index corresponding to sideward direction
        int iOut;           //! In a 3-vector, this is index corresponding to outward direction 
        int iLong;          //! In a 3-vector, this is index corresponding to longitudinal direction 
        double S[3];        //! Unit vector in sideward direction
        double O[3];        //! Unit vector in outward direction
        double L[3];        //! Unit vector in longitudinal direction (set in constructor)
        double q_lab[4];    //! Relative 4-momentum \f$ q = \frac{1}{2}(p_1-p_2) \f$ in lab frame
        double q_cm[4];     //! Relative 4-momentum \f$ q = \frac{1}{2}(p_1-p_2) \f$ in COM frame
        double P_lab[4];    //! Total pair momentum \f$ P = (p_1+p_2) \f$ in lab frame
        double beta[4];     //! Boost 4-velocity from lab to COM frame
        double r_lab[4];    //! Space-time separation, \f$ r=(r_0,\vec{r}) \f$ , of pair emission in lab frame
        double r_cm[4];     //! Space-time separation, \f$ r=(r_0,\vec{r}) \f$ , of pair emission in COM frame
        double sqr_q_cm;    //! \f$ q^2 = q_0^2-\vec{q}^2 = -\vec{q}^2 \f$ in COM frame
        double rinv;        //! \f$ |\vec{r}| \f$ in COM frame
        double qinv;        //! $\f |\vec{q}| \f$ in COM frame
};

#endif
