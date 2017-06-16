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
#ifndef __OSCAR_H__
#define __OSCAR_H__

#include "parametermap.h"
#include "sou1d_histo.h"
#include "sou3d_histo.h"
#include <vector>
#include <string>
#include <istream>

using namespace std;

//! One OSCAR-style particle
struct COSCARLine{ 
    //! particle ID, see the Monte-Carlo numbering scheme in the Particle Data Group Review of Particle Properties
    int pid;        
    //! index of the particle in the OSCAR file (unused)
    int index;      
    //! freeze-out position in fm (or fm/c)
    double x[4];    
    //! freeze-out momentum in GeV/c (of GeV)
    double p[4];    
    //! mass of particle in GeV (unused)
    double mass;    
};

vector<COSCARLine> readOSCARFile(string oscarFile);
vector<COSCARLine> filterPID(const vector<COSCARLine> &lines, int pid);
CSourceFtn1dHisto getOSCARSource1d(vector<COSCARLine> lines, parameterMap p);
CSourceFtn3dHisto getOSCARSource3d(vector<COSCARLine> lines, parameterMap p);
istream& operator>>(istream& i, COSCARLine& line);

#endif
