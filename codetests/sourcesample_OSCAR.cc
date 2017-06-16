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
#include <cstdlib>
#include <cmath>
#include <cstdio>

//
// This code calls the wrong function prototype for CSourceCalc_OSCAR::CalcS ??
//
// DAB 11/5/2009
//

using namespace std;

#include "random.h"
#include "sf.h"
#include "cheezyparser.h"
#include "sourcecalc.h"

int main(){
  int ir,nsample;
  char OSCARfilename[160],sdirname[160];
  char ArrayParsFilename[160];
  CCHArray *Asource;

  // Create Array for storing source info
  sprintf(ArrayParsFilename,"test_data/sourcesample_OSCAR/apars_oscar.dat");
  Asource=new CCHArray(ArrayParsFilename);

  // Initialize Source Calc Object
  CSourceCalc_OSCAR *scalc;
  scalc=new CSourceCalc_OSCAR();
  parameter::set(scalc->spars,"OSCARfilename",
		 "../../gromit_hbt/output_data/bjoscar_pionsonly.output.tmp");
  parameter::set(scalc->spars,"Pt",100);
  scalc->SetIDs(211,211);
  cout << scalc->spars << endl;

  // Calculate source array
  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

  // Write source info to file
  sprintf(sdirname,"sdata/OSCAR_pionsonly_kt%g\0",
	  parameter::getD(scalc->spars,"Pt",-200)/2);
  Asource->WriteAllA(sdirname);
}

