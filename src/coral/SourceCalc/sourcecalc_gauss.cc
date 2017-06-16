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
#ifndef __INCLUDE_SOURCECALC_GAUSS__
#define __INCLUDE_SOURCECALC_GAUSS__
#include "sourcecalc.h"

using namespace std;

CSourceCalc_Gaussian::CSourceCalc_Gaussian(){
  InitSPars();
	randy=new CRandom(-1234);
}

void CSourceCalc_Gaussian::InitSPars(){
  // DEFAULT VALUES
  parameter::set(spars,"lambda",1.0);
  parameter::set(spars,"Rx",4);
  parameter::set(spars,"Ry",4);
  parameter::set(spars,"Rz",4);
  parameter::set(spars,"Xoff",0.0);
  parameter::set(spars,"Yoff",0.0);
  parameter::set(spars,"Zoff",0.0);
  parameter::set(spars,"Euler_Phi",0.0);
  parameter::set(spars,"Euler_Theta",0.0);
  parameter::set(spars,"Euler_Psi",0.0);  
}

void CSourceCalc_Gaussian::SetSPars(double lambdaset,
				    double Rxset,double Ryset,double Rzset,
				    double Xoffset,double Yoffset,
				    double Zoffset,double ephiset,
				    double ethetaset,double epsiset){
  parameter::set(spars,"lambda",lambdaset);
  parameter::set(spars,"Rx",Rxset);
  parameter::set(spars,"Ry",Ryset);
  parameter::set(spars,"Rz",Rzset);
  parameter::set(spars,"Xoff",Xoffset);
  parameter::set(spars,"Yoff",Yoffset);
  parameter::set(spars,"Zoff",Zoffset);
  parameter::set(spars,"Euler_Phi",ephiset);
  parameter::set(spars,"Euler_Theta",ethetaset);
  parameter::set(spars,"Euler_Psi",epsiset);
}

void CSourceCalc_Gaussian::SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,double Xoffset, double Yoffset,double Zoffset){
  parameter::set(spars,"lambda",lambdaset);
  parameter::set(spars,"Rx",Rxset);
  parameter::set(spars,"Ry",Ryset);
  parameter::set(spars,"Rz",Rzset);
  parameter::set(spars,"Xoff",Xoffset);
  parameter::set(spars,"Yoff",Yoffset);
  parameter::set(spars,"Zoff",Zoffset);
}
void CSourceCalc_Gaussian::SetSPars(double lambdaset,
				    double Rxset,double Ryset,double Rzset){
  parameter::set(spars,"lambda",lambdaset);
  parameter::set(spars,"Rx",Rxset);
  parameter::set(spars,"Ry",Ryset);
  parameter::set(spars,"Rz",Rzset);
}

void CSourceCalc_Gaussian::CalcAlpha(double **alpha,CCHArray *A){
  const double PI=4.0*atan(1.0);
  double lambda,Rx,Ry,Rz,Xoff,Yoff,Zoff;
  double phi,theta,psi,ctheta,cphi,cpsi,stheta,sphi,spsi;
  int i,j,k;
  double beta[4][4],gamma[4][4],off[4];

  lambda=parameter::getD(spars,"lambda",1.0);
  Rx=parameter::getD(spars,"Rx",4);
  Ry=parameter::getD(spars,"Ry",4);
  Rz=parameter::getD(spars,"Rz",4);

  Xoff=parameter::getD(spars,"Xoff",0);
  Yoff=parameter::getD(spars,"Yoff",0);
  Zoff=parameter::getD(spars,"Zoff",0);

  phi=parameter::getD(spars,"Euler_Phi",0);
  theta=parameter::getD(spars,"Euler_Theta",0);
  psi=parameter::getD(spars,"Euler_Psi",0);
  
  if(A->GetXSYM() && (fabs(Xoff)>1.0E-8 || fabs(phi)>1.0E-8 
	      || fabs(theta)>1.0E-8 || fabs(psi)>1.0E-8)){
    printf("Xsym true, but Xoff, Euler_phi, Euler_Theta or Euler_Psi !=0\n");
    exit(1);
  }
  if(A->GetYSYM() && (fabs(Yoff)>1.0E-8 || fabs(phi)>1.0E-8 ||
	      fabs(psi)>1.0E-8)){
    printf("Ysym true, but Yoff, Euler_phi, or Euler_Psi !=0\n");
    exit(1);
  }
  if(A->GetZSYM() && fabs(theta)>1.0E-8){
    printf("Zsym true, but Zoff or Euler_Theta !=0\n");
    exit(1);
  }

  ctheta=cos(theta); cphi=cos(phi); cpsi=cos(psi);
  stheta=sin(theta); sphi=sin(phi); spsi=sin(psi);
  double U[4][4],Udagger[4][4];
  U[1][1]=ctheta*cphi*cpsi-sphi*spsi;
  U[1][2]=ctheta*sphi*cpsi+cphi*spsi;
  U[1][3]=-stheta*cpsi;
  U[2][1]=-ctheta*cphi*spsi-sphi*cpsi;
  U[2][2]=-ctheta*sphi*spsi+cphi*cpsi;
  U[2][3]=stheta*spsi;
  U[3][1]=stheta*cphi;
  U[3][2]=stheta*sphi;
  U[3][3]=ctheta;
  
  for(i=0;i<=3;i++){
    for(j=0;j<=3;j++){
      Udagger[i][j]=U[j][i];
      alpha[i][j]=beta[i][j]=gamma[i][j]=0.0;
    }
  }
  
  gamma[1][1]=-0.25/(Rx*Rx);
  gamma[2][2]=-0.25/(Ry*Ry);
  gamma[3][3]=-0.25/(Rz*Rz);

  for(i=1;i<=3;i++){    
    for(j=1;j<=3;j++) beta[i][j]+=gamma[i][i]*Udagger[i][j];
  }
  for(i=1;i<=3;i++){
    for(j=1;j<=3;j++){
      for(k=1;k<=3;k++){
	alpha[i][j]+=U[i][k]*beta[k][j];
      }
    }
  } 
  
  off[1]=Xoff;
  off[2]=Yoff;
  off[3]=Zoff;
  alpha[0][0]=-log(Rx*Ry*Rz*pow(4.0*PI,1.5)/lambda);
  for(i=1;i<=3;i++){
    alpha[0][i]=0.0;
    for(j=1;j<=3;j++){
      alpha[0][i]-=2*alpha[i][j]*off[j];
      alpha[0][0]+=alpha[i][j]*off[i]*off[j];
    }
    alpha[i][0]=alpha[0][i];
  }

}  

void CSourceCalc_Gaussian::CalcS(C3DArray *threed){
	int ix,iy,iz,isx,isy,isz;
	int nsx,nsy,nsz;
	int nxmax=threed->GetNXMAX();
	int nymax=threed->GetNYMAX();
	int nzmax=threed->GetNZMAX();
	double delx=threed->GetDELX();
	double dely=threed->GetDELY();
	double delz=threed->GetDELZ();
	double x,y,z,svalue,snorm;
	double lambda=parameter::getD(spars,"lambda",1.0);
	double Rx=parameter::getD(spars,"Rx",4);
	double Ry=parameter::getD(spars,"Ry",4);
	double Rz=parameter::getD(spars,"Rz",4);
	double Xoff=parameter::getD(spars,"Xoff",0);
	double Yoff=parameter::getD(spars,"Yoff",0);
	double Zoff=parameter::getD(spars,"Zoff",0);
	const double PI=4.0*atan(1.0);
	snorm=lambda*pow(4.0*PI,-1.5)/(Rx*Ry*Rz);  

	nsx=nsy=nsz=1;
	if(threed->GetXSYM()==0) nsx=2;
	if(threed->GetYSYM()==0) nsy=2;
	if(threed->GetZSYM()==0) nsz=2;
	for(isx=0;isx<nsx;isx++){
		for(ix=0;ix<nxmax;ix++){
			x=(0.5+ix)*delx;
			if(isx==2) x=-x;
			x-=Xoff;
			for(isy=0;isy<nsy;isy++){
				for(iy=0;iy<nymax;iy++){
					y=(0.5+iy)*dely;
					if(isy==1) y=-y;
					y-=Yoff;
					for(isz=0;isz<nsz;isz++){
						for(iz=0;iz<nzmax;iz++){
							z=(0.5+iz)*delz;
							if(isz==1) z=-z;
							z-=Zoff;
							svalue=exp(-0.25*((x*x/(Rx*Rx))+(y*y/(Ry*Ry))+(z*z/(Rz*Rz))));
							threed->SetElement(isx,ix,isy,iy,isz,iz,snorm*svalue);
						}
					}
				}
			}
		}
	}
}

void CSourceCalc_Gaussian::CalcS(CCHArray *A){
  int i,ir;
  int GNMAX,NRADIAL;
  double **alpha;
  double r,DELR,norm;
  CCHArray *b,*bb,*C,*oldC;
  bool XSYM,YSYM,ZSYM;
  XSYM=A->GetXSYM();
  YSYM=A->GetYSYM();
  ZSYM=A->GetZSYM();
  GNMAX=40;
  NRADIAL=A->GetNRADIAL();
  DELR=A->GetRADSTEP();
  C=new CCHArray(GNMAX,1,DELR,XSYM,YSYM,ZSYM);
  oldC=new CCHArray(GNMAX,1,DELR,XSYM,YSYM,ZSYM);
  b=new CCHArray(2,1,DELR,XSYM,YSYM,ZSYM);
  bb=new CCHArray(2,1,DELR,XSYM,YSYM,ZSYM);
  alpha=new double *[4];
  for(i=0;i<4;i++) alpha[i]=new double[4];

  CalcAlpha(alpha,A);
  norm=exp(alpha[0][0]);

  for(ir=0;ir<NRADIAL;ir++){
    r=DELR*(0.5+ir);

    b->SetElement(0,0,0,0,0.0);

    if(!XSYM)
      b->SetElement(1,0,0,0,r*alpha[1][0]);
    if(!YSYM)
      b->SetElement(0,1,0,0,r*alpha[2][0]);
    if(!ZSYM)
      b->SetElement(0,0,1,0,r*alpha[3][0]);

    b->SetElement(2,0,0,0,r*r*alpha[1][1]);
    if(!XSYM && !YSYM)
      b->SetElement(1,1,0,0,r*r*alpha[1][2]);
    if(!XSYM && !ZSYM)
      b->SetElement(1,0,1,0,r*r*alpha[1][3]);
    b->SetElement(0,2,0,0,r*r*alpha[2][2]);
    if(!YSYM && !ZSYM)
      b->SetElement(0,1,1,0,r*r*alpha[2][3]);
    b->SetElement(0,0,2,0,r*r*alpha[3][3]);
    //b->PrintArrayFixedIR(0);
    ArrayCalc::Detrace(b,0,bb,0);
    
    ArrayCalc::CalcAExpArrayFromXExpArray(bb,0,A,ir);
    A->ScaleArray(norm,ir);
    A->FillRemainderX(ir);
  }
  for(i=0;i<4;i++) delete [] alpha[i];
  delete [] alpha;
  delete(C);
  delete(oldC);
  delete(b);
  delete(bb);
}

void CSourceCalc_Gaussian::GaussCFCalc(C3DArray *cf3d){
	int ix,iy,iz,isx,isy,isz;
	int nxmax=cf3d->GetNXMAX();
	int nymax=cf3d->GetNYMAX();
	int nzmax=cf3d->GetNZMAX();
	double delqx=cf3d->GetDELX();
	double delqy=cf3d->GetDELY();
	double delqz=cf3d->GetDELZ();
	double qx,qy,qz;
	double lambda=parameter::getD(spars,"lambda",1.0);
	double Rx=parameter::getD(spars,"Rx",4);
	double Ry=parameter::getD(spars,"Ry",4);
	double Rz=parameter::getD(spars,"Rz",4);
	double svalue;  
	isx=isy=isz=0;
	for(ix=0;ix<nxmax;ix++){
		qx=(0.5+ix)*delqx;
		for(iy=0;iy<nymax;iy++){
			qy=(0.5+iy)*delqy;
			for(iz=0;iz<nzmax;iz++){
				qz=(0.5+iz)*delqz;
				svalue=exp((-4*qx*qx*Rx*Rx-4*qy*qy*Ry*Ry-4*qz*qz*Rz*Rz)/(HBARC*HBARC));
				cf3d->SetElement(isx,ix,isy,iy,isz,iz,lambda*svalue);
			}
		}
	}
}

#endif
