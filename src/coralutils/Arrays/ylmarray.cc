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
#ifndef __INCLUDE_YLMARRAY
#define __INCLUDE_YLMARRAY
#include "arrays.h"

// MEMBER FUNCTIONS of CYlmArray

CYlmArray::CYlmArray(int LMAXset,int NRADIALset){
  int L,M,ir;
  LMAX=LMAXset;
  NRADIAL=NRADIALset;
  ylm=new complex<double> **[LMAX+1];
  for(L=0;L<=LMAX;L++){
    ylm[L]=new complex<double> *[L+1];
    for(M=0;M<=L;M++){
      ylm[L][M]=new complex<double> [NRADIAL];
      for(ir=0;ir<NRADIAL;ir++) ylm[L][M][ir]=0.0;
    }
  }
}

CYlmArray::~CYlmArray(){
  int L,m;
  for(L=0;L<=LMAX;L++){
    for(m=0;m<=L;m++) delete [] ylm[L][m];
    delete ylm[L];
  }
  delete ylm;
}

int CYlmArray::GetLMAX(){
  return LMAX;
}

complex<double> CYlmArray::GetElement(int L,int M,int ir){
  return ylm[L][M][ir];
}

void CYlmArray::SetElement(int L,int M,int ir,complex<double> element){
  ylm[L][M][ir]=element;
}

void CYlmArray::IncrementElement(int L,int M,int ir,complex<double> increment){
  ylm[L][M][ir]+=increment;
}

void CYlmArray::ScaleArray(double scalefactor){
  int L,M,ir;
  for(L=0;L<=LMAX;L++){
    for(M=0;M<=L;M++){
      for(ir=0;ir<NRADIAL;ir++) ylm[L][M][ir]*=scalefactor;
    }
  }
}

void CYlmArray::ScaleArray(double scalefactor,int ir){
  int L,M;
  for(L=0;L<=LMAX;L++){
    for(M=0;M<=L;M++){
      ylm[L][M][ir]*=scalefactor;
    }
  }
}

void CYlmArray::ZeroArray(){
  int L,M,ir;
  for(L=0;L<=LMAX;L++){
    for(M=0;M<=L;M++){
      for(ir=0;ir<NRADIAL;ir++) ylm[L][M][ir]=0.0;
    }
  }
}

void CYlmArray::ZeroArray(int ir){
  int L,M;
  for(L=0;L<=LMAX;L++){
    for(M=0;M<=L;M++){
      ylm[L][M][ir]=0.0;
    }
  }
}

void CYlmArray::PrintArrayFixedIR(int ir){
  int L,M;

  printf("\n______________________________________\n");
  printf(" L \\ M :");
  for(M=0;M<=LMAX;M++) printf("         %2d             ",M);
  printf("\n");
  for(L=0;L<=LMAX;L++){
    printf(" %3d : ",L);
    for(M=0;M<=LMAX;M++) printf("(%10.3e,%10.3e) ",
				real(ylm[L][M][ir]),imag(ylm[L][M][ir]));
    printf("\n");
  }
  printf("_________________________________________\n");
}

#endif
