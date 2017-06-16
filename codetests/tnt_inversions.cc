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
#include "message.h"
#include "constants.h"

#include "tnt_array1d.h"
#include "tnt_array1d_utils.h"
#include "tnt_array2d.h"
#include "tnt_array2d_utils.h"
#include "linalg.h"
#include "lsqrinvert.h"

using namespace TNT;
using namespace std;

int main( ) {
    
    // Declare variables, define starting points
    int num_constraints = 0;
    int num_data = 4;
    int num_result = 3;
    Array2D< double > kmtx( num_data, num_result, 1.0 );
    Array1D< double > convec( num_constraints, 0.0 );
    Array2D< double > conmtx( num_constraints, num_result, 0.0 );
    Array1D< double > data( num_data, 1.0 );
    Array2D< double > data_covmtx( num_data, num_data, 1.0 );
    Array1D< double > result( num_result, 0.0 );
    Array2D< double > result_covmtx( num_result, num_result, 0.0 );
    Array1D< double > recons( num_data, 0.0 );
    Array2D< double > recons_covmtx( num_data, num_data, 0.0 );
    
    // Invert
    if ( num_constraints > 0 ) {
        cout << "    Using constraints"<<endl;
        CLSqrInvertSVDLagrange inverter( kmtx, conmtx, convec );
        inverter.solve( data, data_covmtx );
        result        = inverter.model();
        result_covmtx = inverter.covmodel();
        cout << "    Lagrange Multiplier vector:\n";
        cout << inverter.lagrange_multipliers() << endl;
    } else 
        LeastSquaresInvert( data, data_covmtx, kmtx, result, result_covmtx );

    // Restore
    recons = kmtx * result;
    recons_covmtx = matmult( matmult( kmtx, result_covmtx ), transpose( kmtx ) );
    
    // Print out results
    cout << "kmtx:" << kmtx << endl;
    cout << "convec:" << convec << endl;
    cout << "conmtx:" << conmtx << endl;
    cout << "data:" << data << endl;
    cout << "data_covmtx:" << data_covmtx << endl;
    cout << "result:" << result << endl;
    cout << "result_covmtx:" << result_covmtx << endl;
    cout << "recons:" << recons << endl;
    cout << "recons_covmtx:" << recons_covmtx << endl;
    cout << "residual:" << data - recons << endl;
    
}
