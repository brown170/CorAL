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
#include <complex>

#include "arrays.h"
#include "random.h"
#include "misc.h"
#include "parameterMap.h"

#define CHECKMULTDIV
#define CHECKXCONV

using namespace std;

int main(){
  CCHCalc chcalc;
  CRandom random(-1234);
  int lx,ly,lz,m,L,LMAX=6,ir,NRADIAL=1;
  int Lprime,lxprime,lyprime,lzprime;
  double ex,ey,ez,expansion,expA,expAA,expB,expC,ratio,irepeat,lfact,delm,gamma;
  double xarg,expX;
  double *fact;
  bool XSYM,YSYM,ZSYM;
  //XSYM=YSYM=ZSYM=0;
  XSYM=YSYM=ZSYM=0;
  XSYM=1;
  fact=new double[LMAX+1];
  fact[0]=1.0;
  for(m=1;m<=LMAX;m++) fact[m]=double(m)*fact[m-1];

  CCHArray *A,*X,*BB;
  CCHArray *M;
  CYlmArray *Ylm;
  A=new CCHArray(LMAX,NRADIAL,1.0,XSYM,YSYM,ZSYM);
  M=new CCHArray(LMAX,NRADIAL,1.0,XSYM,YSYM,ZSYM);
  Ylm=new CYlmArray(LMAX,NRADIAL);

  for(ir=0;ir<NRADIAL;ir++){
    printf("_______________________________________________________________\n");
    printf("__________________________ ir=%d ______________________________\n",ir);
  TRY_AGAIN:
    ex=1.0-2.0*random.ran();
    ey=1.0-2.0*random.ran();
    if(ex*ex+ey*ey>1.0) goto TRY_AGAIN;
    ez=sqrt(1.0-ex*ex-ey*ey);
    printf("ex=%g, ey=%g, ez=%g, e^2=%g\n",ex,ey,ez,ex*ex+ey*ey+ez*ez);
  
    //printf("Enter lx, ly and lz : ");
    //scanf("%d %d %d",&lx,&ly,&lz);
    lx=2;
    ly=0;
    lz=0;
    L=lx+ly+lz;
    lfact=1.0;
    for(m=1;m<=L;m++) lfact=lfact*(2.0*m+1)/double(m);

    expansion=chcalc.GetMFromE(lx,ly,lz,ex,ey,ez);
    printf("M=%12.5e, ",expansion);
    
    expansion=chcalc.GetAFromE(lx,ly,lz,ex,ey,ez);
    printf("A=%12.5e\n",lfact*expansion);
    
    //

    M->ZeroArray(ir);
    M->IncrementMArrayFromE(ex,ey,ez,1.0,ir);
    printf("M=%12.5e, ",M->GetElement(lx,ly,lz,ir));

    A->ZeroArray(ir);
    A->IncrementAExpArrayFromE(ex,ey,ez,1.0,ir);
    //A->CalcAExpArrayFromE(ex,ey,ez,ir);
    printf("A=%12.5e\n",A->GetElement(lx,ly,lz,ir));
      
    //

    expansion=A->GetMElementFromAExpArray(lx,ly,lz,ir);
    printf("M=%12.5e, ",expansion);
    
    expansion=M->GetAExpElementFromMArray(lx,ly,lz,ir);
    printf("A=%12.5e\n",expansion);
    
    //

    printf("check Moment<->AArray conversion\n");

    for(irepeat=0;irepeat<3;irepeat++){
      ArrayCalc::CalcMArrayFromAExpArray(A,ir,M,ir);
      printf("M=%12.5e, ",M->GetElement(lx,ly,lz,ir));
      ArrayCalc::CalcAExpArrayFromMArray(M,ir,A,ir);
      printf("A=%12.5e\n",A->GetElement(lx,ly,lz,ir));
    }
    expansion=A->GetElement(lx,ly,lz,ir);
    //

    printf("_________ check ylm conversions ___________\n");
    for(irepeat=0;irepeat<3;irepeat++){
      ArrayCalc::CalcYlmExpArrayFromAExpArray(A,ir,Ylm,ir);
      ArrayCalc::CalcAExpArrayFromYlmExpArray(Ylm,ir,A,ir);
      printf("                A=%12.5e\n",A->GetElement(lx,ly,lz,ir));
    }
    printf("__________________\n");
    ratio=A->GetElement(lx,ly,lz,ir)/expansion;
    if(fabs(ratio-1.0)>1.0E-8)
      printf("Expansion off by factor %g (1/factor=%g)\n",ratio,1.0/ratio);

#ifdef CHECKXCONV
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    printf("_________ check XExp conversions ___________\n");
    X=new CCHArray(4,NRADIAL,1.0,XSYM,YSYM,ZSYM);
    BB=new CCHArray(36,1,1.0,XSYM,YSYM,ZSYM);
    X->Randomize(0.4,ir); 
    X->FillRemainderX(ir);
    X->SetElement(0,0,0,ir,1.0);
    printf("                X=%12.5e\n",X->GetElement(lx,ly,lz,ir));
    //X->PrintArrayFixedIR(2,ir);
    for(irepeat=0;irepeat<3;irepeat++){
      BB->ZeroArray();
      ArrayCalc::CalcAExpArrayFromXExpArray(X,ir,BB,0);
      xarg=X->AExpand(ex,ey,ez,ir);
      expX=BB->AExpand(ex,ey,ez,0);
      printf("exp(X)=%g, BB=%g\n",exp(xarg),expX);
      X->ZeroArray();
      ArrayCalc::CalcXExpArrayFromAExpArray(BB,0,X,ir);
      printf("                X=%12.5e\n",X->GetElement(lx,ly,lz,ir));
      //X->PrintArrayFixedIR(2,ir);
    }
    printf("__________________\n");
    //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#endif

    printf("Now check moment expansion\n");
    M->ZeroArray(ir);
    A->ZeroArray(ir);
    M->Randomize(1.0,ir);

    expansion=M->AExpand(ex,ey,ez,ir);
    printf("Expansion yields %g ",expansion);

    ArrayCalc::Detrace(M,ir,A,ir);
    expansion=A->AExpand(ex,ey,ez,ir);
    printf("=? %g\n",expansion);
  
#ifdef CHECKMULTDIV
    printf("_______________________________________\n");
    printf("Now check Multiplication and Division\n");
    CCHArray *B,*C,*AA;


    AA=new CCHArray(2*LMAX,NRADIAL,1.0,XSYM,YSYM,ZSYM);
    B=new CCHArray(LMAX,NRADIAL,1.0,XSYM,YSYM,ZSYM);
    B->Randomize(1.0,ir);
    B->FillRemainderX(ir);
    B->SetElement(0,0,0,ir,100.0);
    C=new CCHArray(LMAX,NRADIAL,1.0,XSYM,YSYM,ZSYM);
    C->Randomize(1.0,ir);
    C->FillRemainderX(ir);
    C->SetElement(0,0,0,ir,1.0);

    expB=B->AExpand(ex,ey,ez,ir);
    expC=C->AExpand(ex,ey,ez,ir);
    ArrayCalc::MultiplyArrays(B,ir,C,ir,AA,ir);
    expAA=AA->AExpand(ex,ey,ez,ir);
    printf("Check Multiplication: expB*expC=%g =? %g=expAA\n",expB*expC,expAA);

    printf("Before Mult and Divide expC=%g\n",expC);
    //C->PrintArrayFixedIR(4,ir);
    ArrayCalc::DivideArrays(AA,ir,B,ir,C,ir);
    expC=C->AExpand(ex,ey,ez,ir);
    printf("After Mult and Divide, expC=%g\n",expC);
    //C->PrintArrayFixedIR(4,ir);
#endif
  }

  return 0;
}
