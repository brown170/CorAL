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
#ifndef __CORAL_CONSTANTS_H__
#define __CORAL_CONSTANTS_H__

#include <cmath>
#include <complex>

const double ZERO         = 0.0;
const double HBARC        = 197.3269602;       // hbar times c
const double ALPHA        = 1.0/137.03599976;   // fine structure constant
const double PI           = 3.141592653589793238462643383279;  
const double SQRTPI       = 1.772453850905516;
const double SQRTFOURPI   = 3.544907702;
const double DEGRAD       = 57.29577951308232;
const double AMU          = 931.494;          // atomic mass unit
const double ProtonMass   = 938.272;
const double KaonMass     = 493.677;
const double PionMass     = 139.57018;
const double Pion0Mass    = 134.9766;
const double LambdaMass   = 1115.7;
const double NeutronMass  = 939.565;
const double RhoMass      = 771.1;
const double XiMass       = 1321.3;
const double XiStarMass   = 1530.0;
const std::complex< double > ci = std::complex< double >(0.0,1.0);

#endif
