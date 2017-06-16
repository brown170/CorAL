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
#include "coral.h"
#include "xgraph.h"
#include <ctime>
using namespace std;

//
// Do we use this?
//
// DAB 11/5/2009
//

int main(){
  CWaveFunction *wf;

  char title[100];
  CXGraph xgraph(400,400,50,50);
  xgraph.setaxes(0.0,-1.0,80.0,1.0);
  xgraph.drawaxes();
  sprintf(title,"PP Correlations, R=3.0 fm\0");
  xgraph.drawtext(title,40.0,0.8);
  xgraph.setcolor("cyan");

  double q,r,ctheta,Rx,Ry,Rz,x,y,z;
  int iq,Nq;
  double *c;
  int imc,NMC=100000;
  CRandom random(-time(NULL));
  char parsfilename[120];
  Rx=Ry=3.0; Rz=3.0;
  double Cfile[40],qfile[40];
  FILE *fptr;
  fptr=fopen("testdata/pptest_R3.dat\0","r");
  for(iq=0;iq<40;iq++){
    fscanf(fptr,"%lf %lf",&qfile[iq],&Cfile[iq]);
    printf("iq=%d\n",iq);
    Cfile[iq]=Cfile[iq]-1.0;
  }
  fclose(fptr);
  xgraph.plotline(qfile,Cfile,40);

  xgraph.setcolor("red");
  //printf("Enter N_MonteCarlo : ");
  //scanf("%d",&NMC);

  sprintf(parsfilename,"parameters/wfparameters.dat\0");
  wf=new CWaveFunction_pp_schrod(parsfilename);
  //wf->PrintCdelta(R,R,R);


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
    xgraph.drawdiamond(q,c[iq]-1.0,0.01);
    printf("%5.2f : %g\n",q,c[iq]);
  }
  Misc::Pause();

}

