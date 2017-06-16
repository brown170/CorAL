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
#include <string>
#include <iostream>

using namespace std;

#include "wavefunction.h"
#include "kernel.h"
#include "parametermap.h"
#include "cheezyparser.h"
#include "random.h"

int main(){
  CWaveFunction *wf;
  double q,r,ctheta,Rx,Ry,Rz,x,y,z;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-256);
  Rx=Ry=5.0; Rz=5.0;

  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);

  string parsfilename = "test_data/wfsample/wfparameters.dat";
  wf=new CWaveFunction_pipluspiplus_sqwell( const_cast< char*>(parsfilename.c_str()) );
  //wf->PrintCdelta(Rx,Ry,Rz);

  Nq=wf->GetNQMAX();
  c=new double[Nq];
  for(iq=0;iq<Nq;iq++){
    q=wf->GetQ(iq);
    c[iq]=0.0;
    for(imc=0;imc<NMC;imc++){
      x=Rx*sqrt(2.0)*random.gauss();
      y=Ry*sqrt(2.0)*random.gauss();
      z=Rz*sqrt(2.0)*random.gauss();
      r=sqrt(x*x+y*y+z*z);
      ctheta=z/r;
      //c[iq]+=wf->GetPsiSquared(q,r,ctheta);
      c[iq]+=wf->CalcPsiSquared(iq,r,ctheta);
    }
    c[iq]=c[iq]/double(NMC);
    printf("%5.2f : %g\n",q,c[iq]);
  }

}

