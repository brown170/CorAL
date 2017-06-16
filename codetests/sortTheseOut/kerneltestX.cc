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
using namespace std;

//
// Do we use this?
//
// DAB 11/5/2009
//

#include "coral.h"
#include "xgraph.h"
#include <ctime>

int main(){
  CWaveFunction *wf;
  CKernel *kernel;
	
  char title[100];
  CXGraph xgraph(400,400,50,50);
  xgraph.setaxes(0.0,-1.0,80.0,1.0);
  xgraph.drawaxes();
  sprintf(title,"PP Correlations, R=3.0 fm\0");
  xgraph.drawtext(title,40.0,0.8);
  xgraph.setcolor("cyan");
	
  double q,r,R=3.0,x,y,z,kvalue;
  int iq,Nq;
  double *c;
  int imc,NMC=1000000;
  CRandom random(-time(NULL));
  char parsfilename[120];
  char kdatadirname[200];
  sprintf(kdatadirname,"kdata/pp");
  double Cfile[40],qfile[40];
  FILE *fptr;
  fptr=fopen("testdata/pptest_R3.dat","r");
  for(iq=0;iq<40;iq++){
    fscanf(fptr,"%lf %lf",&qfile[iq],&Cfile[iq]);
    Cfile[iq]=Cfile[iq]-1.0;
  }
  fclose(fptr);
  xgraph.plotline(qfile,Cfile,40);
  xgraph.setcolor("red");
	
  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);
  NMC=100000;
	
  sprintf(parsfilename,"parameters/wfparameters.dat\0");
	
	kernel=new CKernel(parsfilename);
  //__________________________________________________
  // Either read or calc
  //kernel->ReadData(kdatadirname);
	
  wf=new CWaveFunction_pp_schrod(parsfilename);
  kernel->Calc(wf);
	printf("check\n");
  kernel->WriteData(kdatadirname);
  //__________________________________________________
	
  Nq=kernel->GetNQMAX();
  c=new double[Nq];
  for(iq=0;iq<Nq;iq++){
    q=kernel->GetQ(iq);
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
    xgraph.drawdiamond(q,c[iq]-1.0,0.01);
    printf("c(q=%5.2f) = %g\n",q,c[iq]);
  }
  Misc::Pause();
}

