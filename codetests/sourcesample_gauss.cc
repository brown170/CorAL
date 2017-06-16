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

#include "random.h"
#include "sf.h"
#include "cheezyparser.h"
#include "sourcecalc.h"

int main(){
  int ir;
  double x,y,z,ctheta,phi,guess,r;
  const double PI=4.0*atan(1.0);
  CRandom randy(-12345);
  double Rx=3,Ry=4,Rz=5,Xoff=0.5,Yoff=1.5,Zoff=2.5,lambda=0.6;
  char sdirname[160];
  char ArrayParsFilename[120];
  CCHArray *Asource;
  CSourceCalc_Gaussian *scalc;

  // Create Array to store Cart. Harmonic source info
  sprintf(ArrayParsFilename,"test_data/sourcesample_gauss/apars_gauss.dat");
  Asource=new CCHArray(ArrayParsFilename);
  Asource->PrintPars();

  // Create Source Calc Object
  scalc=new CSourceCalc_Gaussian();
  scalc->SetSPars(lambda,Rx,Ry,Rz,Xoff,Yoff,Zoff);
  cout << scalc->spars << endl;

  // Have Source Calc Object Fill Array
  scalc->CalcS(Asource);
  scalc->NormCheck(Asource);
  scalc->CalcEffGaussPars(Asource);

  // Make 3D Cartesian array from Cart. Harmonic Data
  sprintf(ArrayParsFilename,"test_data/sourcesample_gauss/apars3d.dat");
  C3DArray *threed=new C3DArray(ArrayParsFilename);
  ArrayCalc::Calc3DArrayFromAExpArray(Asource,threed);
  printf("Created 3darray\n");

  // Calc S for some random angles and compare to correct answer
  for(ir=1;ir<Asource->GetNRADIAL();ir+=2){
    ctheta=1.0-2.0*randy.ran();
    phi=2.0*PI*randy.ran();
    r=(0.5+ir)*Asource->GetRADSTEP();
    z=r*ctheta;
    x=r*sqrt(1.0-ctheta*ctheta);
    y=x*sin(phi);
    x=x*cos(phi);
    guess=exp(-0.25*( ((x-Xoff)*(x-Xoff)/(Rx*Rx))
		      +((y-Yoff)*(y-Yoff)/(Ry*Ry))
		      +((z-Zoff)*(z-Zoff)/(Rz*Rz))));
    guess=guess/(Rx*Ry*Rz*pow(4.0*PI,1.5));
    printf("r=%g, S(%g,%g,%g)=%g =? %g =? %g\n",r,x,y,z,
	   Asource->AExpand(x,y,z),guess,
	   threed->GetElement(x,y,z));
  }
  
  // Write array info to file
  sprintf(sdirname,"sample_output/gauss_R%g,%g,%g_off%g,%g,%g",
	  Rx,Ry,Rz,Xoff,Yoff,Zoff);
  Asource->WriteAllA(sdirname);
}

