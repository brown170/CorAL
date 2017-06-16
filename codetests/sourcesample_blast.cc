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

using namespace std;

#include "random.h"
#include "sf.h"
#include "cheezyparser.h"
#include "sourcecalc.h"

int main(){
  int ir,nsample;
  char ArrayParsFilename[120];
  char filename[120];
  CCHArray *Asource,*Xsource;
  sprintf(ArrayParsFilename,"test_data/sourcesample_blast/apars_blast.dat");
  Asource=new CCHArray(ArrayParsFilename);
  sprintf(filename,"Xarray.dat");

  CSourceCalc_Blast *scalc;
  scalc=new CSourceCalc_Blast();
  //printf("Enter nsample (NMC=nsample^2): ");
  //scanf("%d",&nsample);
  nsample=1000;
  parameter::set(scalc->spars,"Nsample",nsample);
  parameter::set(scalc->spars,"etaG",1.0);
  parameter::set(scalc->spars,"Ma",139.59);
  parameter::set(scalc->spars,"Mb",139.59);
  parameter::set(scalc->spars,"DelTau",5.0);
  parameter::set(scalc->spars,"Pt",100.0);
  Asource->PrintPars();
  cout << scalc->spars << endl;

  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

}

