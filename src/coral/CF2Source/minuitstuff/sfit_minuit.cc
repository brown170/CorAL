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
#ifndef __INCLUDE_SFIT_MINUIT_CC__
#define __INCLUDE_SFIT_MINUIT_CC__

#include "sfit_minuit.h"
CMNPars *CCF2S_Minuit::pars=NULL;
double *CCF2S_Minuit::xval=NULL;
double *CCF2S_Minuit::grad=NULL;
C3DArray *CCF2S_Minuit::cexp3D=NULL;
C3DArray *CCF2S_Minuit::cerror3D=NULL;
C3DArray *CCF2S_Minuit::ctheory3D=NULL;
CCHArray *CCF2S_Minuit::source=NULL;
CCHArray *CCF2S_Minuit::ctheory=NULL;
CCHArray *CCF2S_Minuit::cexp=NULL;
CCHArray *CCF2S_Minuit::cerror=NULL;
CKernel *CCF2S_Minuit::kernel=NULL;
CSourceCalc *CCF2S_Minuit::sourcecalc=NULL;
int CCF2S_Minuit::Ncalls=0;
int CCF2S_Minuit::npars=0;
int CCF2S_Minuit::dummy=0;
int CCF2S_Minuit::ndim=0;
int CCF2S_Minuit::lx=0;
int CCF2S_Minuit::ly=0;
int CCF2S_Minuit::lz=0;

void CMNPars::Set(string nameset,double valueset,double errorset,
		  double minset,double maxset){
  strcpy(name,nameset.c_str());
  value=valueset;
  error=errorset;
  min=minset;
  max=maxset;
}

void CCF2S_Minuit::CalcChiSquare(int npar,double* grad,double* fcnval,
				   double* xval,int iflag,void* futil){
  if(ndim==3) CalcChiSquare3D(npar,grad,fcnval,xval,iflag,futil);
  if(ndim==1) CalcChiSquare1D(npar,grad,fcnval,xval,iflag,futil);
}

void CCF2S_Minuit::ViewPars(){
  int ipar,ivarble;
  char name[30];
  double value,error,bnd1,bnd2;
  printf("ipar        name     value       error       bnd1        bnd2    IV(0 for fixed)\n");
  for(ipar=0;ipar<npars;ipar++){
    MNPOUT(ipar+1,name,value,error,bnd1,bnd2,ivarble);
    printf("%2d : %12s %11.4e %11.4e %11.4e %11.4e %d\n",
	   ipar,name,value,error,bnd1,bnd2,ivarble);
  }
}

void CCF2S_Minuit::Mnstat(){
  double Fmin,delFestimate,ErrDef=1.0;
  int npars_variable,npars_total,istat;
  char message[40];
  MNSTAT(Fmin,delFestimate,ErrDef,npars_variable,npars_total,istat);
  printf("Min Val. thus far=%g, Estimate of delF = %g\n",Fmin,delFestimate);
  printf("npars_variable=%d, npars_total=%d\n",npars_variable,npars_total);
  if(istat==0) sprintf(message,"Not calculated at all\n");
  else if(istat==1)
    sprintf(message,"Diagonal approximation only, not accurate\n");
  else if(istat==2)
    sprintf(message,"Full matrix, but forced positive definite\n");
  else if(istat==3)
    sprintf(message,"Full accurate covariance matrix, as good as it gets\n");
  printf("Convergence status = %d (%s)\n",message);
}

void CCF2S_Minuit::SetPar(int ipar,char *name,double value,double error,double min,double max){
  MNPARM(ipar+1,name,value,error,min,max,dummy);
}

void CCF2S_Minuit::FixPar(int ipar){
  char command[MAXLINE];
  int ierr=0;
  //snprintf(command,MAXLINE,"FIX %d\0",ipar+1);
  //MNCOMD(CalcChiSquare,command,dummy,0);
  MNFIXP(ipar+1,ierr);
}

void CCF2S_Minuit::FreePar(int ipar){
  char command[MAXLINE];
  int ierr;
  //snprintf (command, MAXLINE, "FREE %d\0",ipar+1);
  //MNCOMD(CalcChiSquare,command,dummy,0);
  MNFREE(ipar+1);
}

void CCF2S_Minuit::InitMinuit(){
  int i;
  MNINIT(5,6,7);
  for (i=0;i<npars;i++)
    SetPar(i,pars[i].name,pars[i].value,pars[i].error,
	   pars[i].min,pars[i].max);
  printf("Minuit Initialized\n");
}

void CCF2S_Minuit::Scan(int ipar,int npts,double start,double end){
  char command[MAXLINE];
  snprintf(command,MAXLINE,"SCA %d %d %g %g\n",ipar,npts,start,end);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Migrad(){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MIGRAD\0");
  MNCOMD(CalcChiSquare,command,dummy,0);
}
void CCF2S_Minuit::Migrad(int maxcalls){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MIGRAD %d\0",maxcalls);
  MNCOMD(CalcChiSquare,command,dummy,0);
}
void CCF2S_Minuit::Migrad(int maxcalls,double tolerance){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MIGRAD %d %g\0",maxcalls,tolerance);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Minimize(){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MINI\0");
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Minimize(int maxcalls){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MINI %d\0",maxcalls);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Minimize(int maxcalls,double tolerance){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MINI %d %g\0",maxcalls,tolerance);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Simplex(){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "SIM\0");
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Simplex(int maxcalls){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "SIM %d\0",maxcalls);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Simplex(int maxcalls,double tolerance){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "SIM %d %g\0",maxcalls,tolerance);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Minos(){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MINOS\0");
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Minos(int ipar){
  char command[MAXLINE];
  snprintf (command, MAXLINE, "MINOS 100 %d\0",ipar+1);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::StratLevel(int level){
  char command[MAXLINE];
  if(level!=0 && level!=1 && level!=2){
    printf("Attempting to set strategy level to %d - only 0,1,2 allowed\n",
	   level);
    printf("Will set level to 1\n"); 
    level=1;
  }
  snprintf (command, MAXLINE, "SET STR %d\0",level);
  MNCOMD(CalcChiSquare,command,dummy,0);
  
}

void CCF2S_Minuit::SetError(double error){
  char command[MAXLINE];
  snprintf(command,MAXLINE,"SET ERR %g\n",error);
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::ErrorMatrix(){
  char command[MAXLINE];
  snprintf(command,MAXLINE,"SHO COV\n");
  MNCOMD(CalcChiSquare,command,dummy,0);
}

void CCF2S_Minuit::Contour(int iparx,int ipary,int npts,
		      double *xcontour,double *ycontour){
  int nfound_contour;
  MNCONT(CalcChiSquare,iparx+1,ipary+1,npts,xcontour[0],ycontour[0],
	 nfound_contour,0);
}

CCF2S_Minuit::CCF2S_Minuit(){
}


void fcn(int npar, double* grad, double* fcnval,
	 double* xval,int iflag, void* futil){
  printf("why am I here?\n");
  exit(1);
}

void CCF2S_Minuit::CalcChiSquare3D(int npar, double* grad, double* fcnval,
			    double* xval,int iflag, void* futil){
  int ipar;
  Ncalls+=1;
  printf("In CalcChiSquare3D, xval[] = ");
  for(ipar=0;ipar<npars;ipar++) printf("%g,",xval[ipar]);
  printf("\n");
  for(ipar=0;ipar<npars;ipar++){
    parameter::set(sourcecalc->spars,pars[ipar].name,xval[ipar]);
  }
  sourcecalc->CalcS(source);
  S2CF::s2c(source,kernel,ctheory);
  ArrayCalc::Calc3DArrayFromAExpArray(ctheory,ctheory3D);
  *fcnval=CFCalc::GetChiSquared(cexp3D,cerror3D,ctheory3D);
  for(ipar=0;ipar<npars;ipar++) pars[ipar].value=xval[ipar];
  printf("fcnval=%g   (Ncalls=%d)\n",*fcnval,Ncalls);
}

void CCF2S_Minuit::CalcChiSquare1D(int npar, double* grad, double* fcnval,
			      double* xval,int iflag, void* futil){
  int ipar;

  Ncalls+=1;
  printf("In CalcChiSquare1D, xval[] = ");
  for(ipar=0;ipar<npars;ipar++) printf("%g,",xval[ipar]);
  printf("\n");
  for(ipar=0;ipar<npars;ipar++){
    parameter::set(sourcecalc->spars,pars[ipar].name,xval[ipar]);
  }
  sourcecalc->CalcS(lx,ly,lz,source);
  S2CF::s2c(source,kernel,ctheory);
  *fcnval=CFCalc::GetChiSquared(lx,ly,lz,cexp,cerror,ctheory);
  for(ipar=0;ipar<npars;ipar++) pars[ipar].value=xval[ipar];
  printf("fcnval=%g   (Ncalls=%d)\n",*fcnval,Ncalls);
}

#endif

