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
#include <iostream>

using namespace std;

#include "wavefunction.h"
#include "kernel.h"
#include "parametermap.h"
#include "cheezyparser.h"
#include "random.h"


int main(){
  CWaveFunction *wf;
  CKernel *kernel;
  double q,r,R=5.0,x,y,z,kvalue;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-123);
  char parsfilename[120];
  char kdatadirname[200];
  sprintf(kdatadirname,"kdata/pp");

  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);
  NMC=10000;

  sprintf(parsfilename,"test_data/kernelsample/wfparameters.dat");
  wf=new CWaveFunction_pipluspiplus_sqwell(parsfilename);
  kernel=new CKernel(parsfilename);
  //parameterMap pars;  
  //ReadParsFromFile( pars, parsfilename );
  //kernel->Read(pars);
  // Either read or calc
  //kernel->ReadData(kdatadirname);
  kernel->Calc(wf);
  //kernel->Calc_ClassCoul(938.28,493.677,1);

  //kernel->Print();
  kernel->WriteData(kdatadirname);
    printf("kernel read\n");
  Nq=kernel->GetNQMAX();
    printf("Nq=%d\n",Nq);
  c=new double[Nq];
  for(iq=0;iq<Nq;iq++){
    q=(iq+0.5)*kernel->GetDELQ();
    c[iq]=0.0;
    for(imc=0;imc<NMC;imc++){
      x=R*sqrt(2.0)*random.gauss();
      y=R*sqrt(2.0)*random.gauss();
      z=R*sqrt(2.0)*random.gauss();
      r=sqrt(x*x+y*y+z*z);
      kvalue=kernel->GetValue(0,q,r);
      c[iq]+=1.0+kvalue;
      if(kvalue<-1.0 && r>1.0){
	printf("WARNING: r=%g, kernel BELOW ZERO, K=%g!!!!!!!!\n",r,kvalue);
	// exit(1);
      }
    }
    c[iq]=c[iq]/double(NMC);
    printf("c(q=%5.2f) = %g\n",q,c[iq]);
  }
/*
  int ir;
  double ctheta;
  for(iq=0;iq<Nq;iq++){
    q=wf->GetQ(iq);
    for(ir=0;ir<kernel->GetNRMAX();ir++){
      ctheta=1.0-2.0*random.ran();
      r=(0.5+ir)*kernel->GetDELR();
      //printf("iq=%d ir=%d ctheta=%g, ",iq,ir,ctheta);
      printf("phi^2=%g =? %g\n",
	     kernel->GetPsiSquared(q,r,ctheta),wf->CalcPsiSquared(iq,r,ctheta));
    }
  }
*/
}

