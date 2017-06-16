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
#include <sstream>
#include <iomanip>
#include "message.h"
#include "utilities.h"
#include "constants.h"
#include "cheezyparser.h"
#include "kernel.h"

using namespace std;

string getKernelFilename( string datadir, int ell, double q )
{
    stringstream filename;
    filename << datadir 
        << "/ell" << ell 
        << "_q" << setfill('0') << fixed << setw(7) << setprecision(2) << q << ".tmp";
    return filename.str();
}

/**
  *  \brief Constructor for CKernel
  *  \param kparsfilename name of the file that has the parameters to set up the kernel
  */
CKernel::CKernel(string kparsfilename){
  int iq,ir,ell,dell;
  ParsInit(kparsfilename);
  // dell is the difference in angular kernelum
  dell=1;
  if(IDENTICAL) dell=2;
	
  kernel=new double **[ellmax+1];
  P=new double [ellmax+1];
  for(ell=0;ell<=ellmax;ell+=1){
    kernel[ell]=NULL;
    if(ell%dell==0){
      kernel[ell]=new double *[nqmax];
      for(iq=0;iq<nqmax;iq++){
				kernel[ell][iq]=new double[nrmax];
				for(ir=0;ir<nrmax;ir++) kernel[ell][iq][ir]=0.0;
      }
    }
  }
}

/**
*  \brief Destructor for CKernel
*/
CKernel::~CKernel(){
  int ell,iq,dell=1;
	if(IDENTICAL) dell=2;
  delete [] P;
  for(ell=0;ell<=ellmax;ell+=dell){
    for(iq=0;iq<nqmax;iq++) delete [] kernel[ell][iq];
    if(kernel[ell]!=NULL) delete [] kernel[ell];
  }
  delete [] kernel;
  printf("kernel object deleted\n");
  
}

//! Getter method for ellmax
//! \return CKernel::ellmax
int CKernel::GetLMAX(){
  return ellmax;
}

//! Getter method for delr
//! \return CKernel::delr
double CKernel::GetDELR(){
  return delr;
}
double CKernel::GetQ(int iq){
  return (iq+0.5)*delq;
}

//! Getter method for delq
//! \return CKernel::delq
double CKernel::GetDELQ(){
  return delq;
}

//! Getter method for nrmax
//! \return CKernel::nrmax
int CKernel::GetNRMAX(){
  return nrmax;
}

//! Getter method for nqmax
//! \return CKernel::nqmax
int CKernel::GetNQMAX(){
  return nqmax;
}

//! Getter method for IDENTICAL
//! \return CKernel::IDENTICAL
bool CKernel::GetIDENTICAL(){
  return IDENTICAL;
}

/**
*  Performs the linear interpolation to get the value of the kernel at a particular point
*  \param  ell The Legendre order of the requested kernel
*  \param  q   The value of relative momentum at which you want to evaluate the kernel (in MeV/c)
*  \param  r   The value of separation at which you want to evaluate the kernel (in fm)
*  \return Value of the kernel at the requested point
*/
double CKernel::GetValue(int ell,double q,double r) const{
  int iqlow,iqhigh,irlow,irhigh,iq,ir;
  double wrlow,wrhigh,wqlow,wqhigh,answer;
  double qlow,qhigh,rlow,rhigh,qi,ri;
  
  // Perform simple linear extrapolation
  if((IDENTICAL && pow(-1.0,ell)<0) || q>delq*nqmax || r>delr*nrmax) answer=0.0;
  else{
    iq=lrint(floor(q/delq));
    qi=(0.5+iq)*delq;
    if(q<qi){
      iqlow=iq-1;
      if(iqlow<0) iqlow=0;
      iqhigh=iq;
      qhigh=(0.5+iq)*delq;
      wqhigh=fabs(delq-fabs(qhigh-q))/delq;
      if(wqhigh>1.0) wqhigh=1.0;
      wqlow=1.0-wqhigh;
    }
    else{
      iqlow=iq;
      iqhigh=iq+1;
      if(iqhigh>=nqmax) iqhigh=iqlow;
      qlow=(iq+0.5)*delq;
      wqlow=fabs(delq-fabs(q-qlow))/delq;
      if(wqlow>1.0) wqlow=1.0;
      wqhigh=1.0-wqlow;
    }
		
    ir=int(floor(r/delr));
    ri=(0.5+ir)*delr;
    if(r<ri){
      irlow=ir-1;
      if(irlow<0) irlow=0;
      irhigh=ir;
      rhigh=(0.5+ir)*delr;
      wrhigh=fabs(delr-fabs(rhigh-r))/delr;
      //wrhigh=1.0;
      if(wrhigh>1.0) wrhigh=1.0;
      wrlow=1.0-wrhigh;
    }
    else{
      irlow=ir;
      irhigh=ir+1;
      if(irhigh>=nrmax) irhigh=irlow;
      rlow=(ir+0.5)*delr;
      wrlow=fabs(delr-fabs(r-rlow))/delr;
      if(wrlow>1.0) wrlow=1.0;
      wrhigh=1.0-wrlow;
    }
		
    if(fabs(wrlow+wrhigh-1.0)>1.0E-8 || fabs(wqlow+wqhigh-1.0)>1.0E-6){
      printf("weights do not add up in GetValue(ell,q,r)\n");
      exit(1);
    }
		
    answer=wqlow*wrlow*GetValue(ell,iqlow,irlow)
		+wqlow*wrhigh*GetValue(ell,iqlow,irhigh)
		+wqhigh*wrlow*GetValue(ell,iqhigh,irlow)
		+wqhigh*wrhigh*GetValue(ell,iqhigh,irhigh);
		
    if(ell==0 && r>1.0 && answer<-1.0){
			printf("__________________________________________________\n");
      printf("answer below -1, =%g, q=%g, r=%g, iq=%d, ir=%d, nrmax=%d\n",
				answer,q,r,iq,ir,nrmax);
			printf("wqlow=%g, wqhigh=%g, wrlow=%g, wrhigh=%g\n",wqlow,wqhigh,wrlow,wrhigh);
      //exit(1);
			printf("__________________________________________________\n");
		}
  }
	return answer;
}

/**
* \brief  The value of a particular point in the interpolation table
* \param  ell The Legendre order of the table in question
* \param  iq  The q grid index
* \param  ir  The r grid index
* \return The value of the interpolation table at the specified index on the grid
*/
double CKernel::GetValue(int ell,int iq,int ir) const{
  if(iq>=nqmax||ir>=nrmax || iq<0 || ir<0) return 0.0;
  else{
    /* 
		if(ell==0 && ir>=2&& kernel[ell][iq][ir]<-1.01){
		printf("Kernel Warning: ell=0, iq=%d, ir=%d, kernel=%g\n",iq,ir,
		kernel[ell][iq][ir]);
		//exit(1);
		} */
    return kernel[ell][iq][ir];
  }
}

/**
* \brief Read in the parameters from a parameterMap
* \param parameters The parameterMap containing the parameters (see below)
* \return True==Success, False==Failure
*
* The parameters in the parameterMap that this code uses are as follows:
*   - \c nqmax     The number of points in the q direction in the kernel interpolation table
*   - \c nrmax     The number of points in the r direction in the kernel interpolation table
*   - \c kellmax   The maximum order of the Legendre polynomial expansion of the kernel to compute
*   - \c delq      The step size in the q direction in the kernel interpolation table
*   - \c delr      The step size in the r direction in the kernel interpolation table
*   - \c IDENTICAL Bool to denote whether this kernel is for identical pairs
*/
bool CKernel::Read(const parameterMap& parameters){
	nqmax=parameter::getI(parameters,"NQMAX",200);
	nrmax=parameter::getI(parameters,"NRMAX",500);
	ellmax=parameter::getI(parameters,"KLMAX",4);
	delq=parameter::getD(parameters,"DELQ",0.5);
	delr=parameter::getD(parameters,"DELR",0.1);
	IDENTICAL=parameter::getB(parameters,"IDENTICAL",true);
	return true;
}


/**
* \brief Writes the parameters to a parameterMap
* \param parameters The parameterMap containing the parameters (see below)
* \return True==Success, False==Failure
*
* The parameters in the parameterMap that this code uses are as follows:
*   - \c nqmax     The number of points in the q direction in the kernel interpolation table
*   - \c nrmax     The number of points in the r direction in the kernel interpolation table
*   - \c kellmax   The maximum order of the Legendre polynomial expansion of the kernel to compute
*   - \c delq      The step size in the q direction in the kernel interpolation table
*   - \c delr      The step size in the r direction in the kernel interpolation table
*   - \c IDENTICAL Bool to denote whether this kernel is for identical pairs
*/
bool CKernel::Write( parameterMap& parameters){
	parameter::set(parameters,"NQMAX",nqmax);
	parameter::set(parameters,"NRMAX",nrmax);
	parameter::set(parameters,"KLMAX",ellmax);
	parameter::set(parameters,"DELQ",delq);
	parameter::set(parameters,"DELR",delr);
	parameter::set(parameters,"IDENTICAL",IDENTICAL);
	return true;
}


void CKernel::ParsInit(string kparsfilename){
  parameterMap parameters;
  if ( kparsfilename!=string("") ) {
    try {
			parameter::ReadParsFromFile(parameters,kparsfilename);
    } catch ( const CMessage& ) {}
  }
  Read(parameters);
	
  printf("    Parameters for kernel: \n");
  printf("        delq set to %g\n",delq);
  printf("        nqmax set to %d\n",nqmax);
  printf("        delr set to %g\n",delr);
  printf("        nrmax set to %d\n",nrmax);
  printf("        KLMAX=%d\n",ellmax);
  printf("        IDENTICAL set to %d\n",IDENTICAL);
//  printf("__________________________________________\n");
}

/**
* \brief Read in the interpolation table from a directory
* \param datadir The directory containing the interpolation table.  
*
* The data in datadir is stored in one file per q.  Each file contains 
* an interpolation table in r for that value of q.
*/
void CKernel::ReadData(string datadir){
  double q,delr0;//,r;
  int iq,ir,ell,nrmax0,dell;
  FILE *infile;
  dell=1;
  string filename;
  if(IDENTICAL) dell=2;
  // Check if all the needed files are there
  for(ell=0;ell<=ellmax;ell+=dell){
    bool this_l_OK = true;
    for(iq=0;iq<nqmax;iq++){
      q=(0.5+iq)*delq;
      filename = getKernelFilename( datadir, ell, q );
      if (! file_exists( filename ) ) this_l_OK = false;
    }
    if ( !this_l_OK )
    {
      MESSAGE << "Missing files for l = "<< ell<<", truncating kernel at lmax = "<<ell-dell<< ENDM_WARN;
      ellmax = ell-dell;
      break;
    }
  } 
  // Load the files
  for(ell=0;ell<=ellmax;ell+=dell){
    for(iq=0;iq<nqmax;iq++){
      q=(0.5+iq)*delq;
      filename = getKernelFilename( datadir, ell, q );
      infile=fopen(filename.c_str(),"r");
      if (infile==NULL) throw MESSAGE<<"Opening kernel data file "<<filename<<" failed!"<<ENDM_SEVERE;
      fscanf(infile,"%d %lf",&nrmax0,&delr0);
      if(nrmax0<nrmax || fabs(delr-delr0)>1.0E-8){
				printf("Inconsistent values for nrmax or delr in data files!\n");
				exit(1);
      }
      for(ir=0;ir<nrmax;ir++){
				fscanf(infile,"%lf",&kernel[ell][iq][ir]);
      }
      fclose(infile);
    }
  }
}

/** 
* \brief Writes the interpolation table to bunch of files in a directory
* \param datadir The directory containing the interpolation table.  
*
* The data in datadir is stored in one file per q.  Each file contains 
* an interpolation table in r for that value of q.
*/
void CKernel::WriteData(string datadir){
  double q,r;
  int ell,iq,ir,dell;
  char filename[100];
  FILE *outfile;
  char shellcommand[120];
  //sprintf(shellcommand,"#!sh; if [ ! -e  %s ]; then mkdir -p %s; fi; \0",
  //  datadir,datadir);
  sprintf(shellcommand,"mkdir -p %s",datadir.c_str());
	
  system(shellcommand);
	
  dell=1;
  if(IDENTICAL) dell=2;
  for(ell=0;ell<=ellmax;ell+=dell){
    for(iq=0;iq<nqmax;iq++){
      q=(0.5+iq)*delq;
      sprintf(filename,"%s/ell%d_q%07.2f.tmp",datadir.c_str(),ell,q);
      printf("For q=%g, Will write to %s\n",q,filename);
      outfile=fopen(filename,"w");
      fprintf(outfile,"%d %g\n",nrmax,delr);
      for(ir=0;ir<nrmax;ir++){
				r=(0.5+ir)*delr;
				fprintf(outfile,"%17.10e\n",kernel[ell][iq][ir]);
      }
      fclose(outfile);
    }
  }
}

/** 
* \brief Simple routine to print out the entire interpolation table to stdout for debugging
*/
void CKernel::Print(){
  double q,r;
  int ell,iq,ir,dell;//,nrmax0;
	
  dell=1;
  if(IDENTICAL) dell=2;
  printf("ellmax=%d, nqmax=%d, nrmax=%d\n",ellmax,nqmax,nrmax);
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    printf("______________ q=%g MeV/c____________________\n",q);
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      printf("%6.2f ",r);
      for(ell=0;ell<=ellmax;ell+=dell){
				printf("%10.3e ",kernel[ell][iq][ir]);
      }
      printf("\n");
    }
  }
}

/**
* \brief Builds the interpolation table, assuming classical Coulomb correlations.  
* \param ma       Mass of particle a, in MeV
* \param mb       Mass of particle b, in MeV
* \param zazb     Charge of particle a times charge or particle b.  Both in units of \f$e\f$.
*
*  THIS DOCUMENTATION NEEDS TO BE FLUSHED OUT!!!!!!
*/
void CKernel::Calc_ClassCoul(double ma,double mb,int zazb){
  double q,r,ctheta;//,delctheta;
  int ell,iq,ir,nu=512;
  double x,ea,eb,u,umax,delu;//,cweight;
	
  for(ell=0;ell<=ellmax;ell++){
    for(iq=0;iq<nqmax;iq++){
      for(ir=0;ir<nrmax;ir++) kernel[ell][iq][ir]=0.0;
    }
  }
  
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      ea=sqrt(q*q+ma*ma);
      eb=sqrt(q*q+mb*mb);
      x=2.0*ea*eb*zazb*197.327/(q*q*r*137.036*(ea+eb));
			
      if(x<1.0){
        umax=2.0*sqrt(1.0-x);
        delu=umax/double(nu);
        for(u=0.5*delu;u<umax;u+=delu){
            ctheta=-1.0+x+sqrt(u*u+x*x);
            for(ell=0;ell<=ellmax;ell++)
            kernel[ell][iq][ir]+=0.5*delu*SpherHarmonics::legendre(ell,ctheta);
        }
        printf("q=%g, r=%g, kernel0=%g =? %g\n",
            q,r,kernel[0][iq][ir],sqrt(1.0-x));
      }
      kernel[0][iq][ir]-=1.0;
    }
  }
}

/**
* \brief Builds the interpolation table, assuming HBT only.  
* 
*  For a pure HBT kernel, the pair relative wavefunction is simply 
*  \f$ \Psi = \frac{1}{\sqrt{2}}\left(e^{iq\cdot r}+e^{-iq\cdot r}\right) \f$
*  Where \f$q = \frac{1}{2}(q_1-q_2)\f$ is the relative four-momentum of the pair and \f$r\f$ is 
*  the space-time separation of the pair at freeze-out.  It is straight-forward to show 
*  that the kernel is
*  \f[
*      K_\ell(r,q) = (-1)^{\ell/2} j_\ell(2qr/\hbar c)
*  \f]
*  This is for even \f$\ell\f$ only.  For odd ones, this is zero.  Here \f$q=|\vec{q}|\f$ 
*  and \f$r=|\vec{r}|\f$, both in the pair rest frame.
*
*/
void CKernel::Calc_PureHBT(){
  int ir,iq,ell,sign;
  double qr,r,q;
	
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      qr=q*r/197.3269602;
      sign=-1;
      for(ell=0;ell<=ellmax;ell+=2){
        sign=-sign;
        kernel[ell][iq][ir]=sign*Bessel::jn(ell,2.0*qr);
//        kernel[ell][iq][ir]=2.0*sign*Bessel::jn(ell,qr);
				//printf("ell=%d ellmax=%d, q=%g, r=%g kernel=%g\n",
				//ell,ellmax,q,r,kernel[ell][iq][ir]);
      }
      kernel[0][iq][ir]-=1.0;
    }
  }
}

/**
* \brief Builds the interpolation table, using the actual CWaveFunction in CKernel::wf
* \param wf pointer to the wavefunction to use to build the kernel
* 
* The table is built up on the regular grid defined by CKernel::nqmax, CKernel::delq,
* CKernel::nrmax, CKernel::delr, by a simple Simpson rule integration over the angle between
* \f$\vec{q}\f$ and \f$\vec{r}\f$:
* \f[
*     K_\ell(q,r) = \frac{1}{2}\int_{-1}^{1} d\mu \left[|\Psi(q,r,\mu)|^2 - 1\right] P_\ell(\mu)
* \f]
* Where \f$\mu = cos(\vec{q}\cdot\vec{r}/qr)\f$.
*/
void CKernel::Calc(CWaveFunction *wf){
  double q,r,ctheta,wf2;
  int ell,dell,iq,ir,ictheta,nctheta=720;
	if(IDENTICAL!=wf->GetIDENTICAL()){
		printf("Creating Kernel with different symmetry (value of IDENTICAL) than wave function\n");
		printf("You probably need to edit parameters file\n");
		exit(1);
	}
	
  dell=1;
  if(IDENTICAL) dell=2;
  for(iq=0;iq<nqmax;iq++){
    q=(0.5+iq)*delq;
    for(ir=0;ir<nrmax;ir++){
      r=(0.5+ir)*delr;
      for(ell=0;ell<=ellmax;ell+=dell)
				kernel[ell][iq][ir]=0.0;
      for(ictheta=0;ictheta<nctheta;ictheta++){
				ctheta=-1.0+2.0*(0.5+ictheta)/double(nctheta);
				wf2=wf->GetPsiSquared(q,r,ctheta);
				for(ell=0;ell<=ellmax;ell+=dell){
					kernel[ell][iq][ir]+=(wf2-1.0)*SpherHarmonics::legendre(ell,ctheta);
				}
      }
      for(ell=0;ell<=ellmax;ell+=dell)
				kernel[ell][iq][ir]=kernel[ell][iq][ir]/double(nctheta);
    }
  }
}

double CKernel::CalcPsiSquared(int iq,int ir,double ctheta){
  int ell,dell=1,nsmall=0;
  double answer=1.0,dela;
  if(IDENTICAL) dell=2;
  ell=0;
  if(ir<nrmax && iq<nqmax){
    while(ell<=ellmax && nsmall<5){
      dela=kernel[ell][iq][ir]*(2.0*ell+1)*P[ell];
      //*SpherHarmonics::legendre(ell,ctheta);
      answer+=dela;
      if(fabs(dela)<1.0E-4) nsmall+=1;
      else(nsmall=0);
      ell+=dell;
    }
  }
  return answer;
}

double CKernel::GetPsiSquared(int iq,int ir,double ctheta){
  double answer=1.0;
  if(ir<nrmax && iq<nqmax){
    CalcP(ctheta);
    answer=CalcPsiSquared(iq,ir,ctheta);
  }
  return answer;
}

double CKernel::CalcPsiSquared(int iq,double r,double ctheta){
  double wa,wb;
  int ira,irb;
  ira=lrint(floor((r/delr)-0.5));
  if(ira<0) ira=0;
  irb=ira+1;
  wb=(r-(ira+0.5)*delr)/delr;
  wa=((irb+0.5)*delr-r)/delr;
  return wa*CalcPsiSquared(iq,ira,ctheta)+wb*CalcPsiSquared(iq,irb,ctheta);
}


double CKernel::GetPsiSquared(int iq,double r,double ctheta){
  double answer=1.0;
  if(iq<nqmax){
    CalcP(ctheta);
    answer=CalcPsiSquared(iq,r,ctheta);
  }
  return answer;
}



double CKernel::GetPsiSquared(double q,double r,double ctheta){
  double wa,wb,answer=1.0;
  int iqa,iqb;
  iqa=lrint(floor((q/delq)-0.5));
  if(iqa<0) iqa=0;
  iqb=iqa+1;
  if(iqb>=nqmax){
    iqb=nqmax-1;
    iqa=iqb-1;
  }
  CalcP(ctheta);
  wb=(q-(iqa+0.5)*delq)/delq;
  wa=((iqb+0.5)*delq-q)/delq;
  answer=wa*CalcPsiSquared(iqa,r,ctheta)+wb*CalcPsiSquared(iqb,r,ctheta);
	
  return answer;
}

void CKernel::CalcP(double ctheta){
  int ell;
  P[0]=1.0;
  P[1]=ctheta;
  for(ell=1;ell<ellmax;ell++){
    P[ell+1]=((2*ell+1)*ctheta*P[ell]-ell*P[ell-1])/double(ell+1);
  }
}    

// Routines Below for CKernelWF

CKernelWF::CKernelWF(string kparsfilename){
  int iq,ir,ictheta;
  ParsInit(kparsfilename);
  // dell is the difference in angular kernelum
  delctheta=2.0/double(nctheta);
  if(IDENTICAL) delctheta=delctheta*0.5;
  
  wfarray=new double **[nqmax];
  for(iq=0;iq<nqmax;iq+=1){
    wfarray[iq]=new double *[nrmax];
    for(ir=0;ir<nrmax;ir++){
      wfarray[iq][ir]=new double[nctheta+1];
      for(ictheta=0;ictheta<=nctheta;ictheta++) wfarray[iq][ir][ictheta]=0.0;
    }
  }
}

void CKernelWF::ParsInit(string kparsfilename){
  parameterMap parameters;
	parameter::ReadParsFromFile(parameters,kparsfilename);
  if ( kparsfilename!="" ) {
    try {
			parameter::ReadParsFromFile(parameters,kparsfilename);
      printf("Read from %s:\n",kparsfilename.c_str());
    } catch ( const CMessage& ) {}
  }
	
  nqmax=parameter::getI(parameters,"NQMAX",25);
  nrmax=parameter::getI(parameters,"NRMAX",120);
  nctheta=parameter::getI(parameters,"NCTHETA",120);
  delq=parameter::getD(parameters,"DELQ",4.0);
  delr=parameter::getD(parameters,"DELR",0.5);
  IDENTICAL=parameter::getB(parameters,"IDENTICAL",0);
	
  printf("reading from %s\n",kparsfilename.c_str());
  printf("  _________ PARAMETERS FOR KERNELWF ________\n");
  printf("delq set to %g\n",delq);
  printf("nqmax set to %d\n",nqmax);
  printf("delr set to %g\n",delr);
  printf("nrmax set to %d\n",nrmax);
  printf("nctheta set to %d\n",nctheta);
  printf("IDENTICAL set to %d\n",IDENTICAL);
  printf("__________________________________________\n");
}


CKernelWF::~CKernelWF(){
  int iq,ir;
  for(iq=0;iq<nqmax;iq++){
    for(ir=0;ir<nrmax;ir++) delete [] wfarray[iq][ir];
    delete [] wfarray[iq];
  }
  delete [] wfarray;
  printf("KernelWF object deleted\n");
  
}

void CKernelWF::Calc(CWaveFunction *wf){
  int iq,ir,ictheta;
  double q,r,ctheta,ps2;
  for(iq=0;iq<nqmax;iq++){
    q=(iq+0.5)*delq;
    for(ir=0;ir<nrmax;ir++){
      r=(ir+0.5)*delr;
      for(ictheta=0;ictheta<=nctheta;ictheta++){
				ctheta=1.0-ictheta*delctheta;
				if(ctheta>1.0) ctheta=1.0;
				ps2=wf->GetPsiSquared(q,r,ctheta);
				if(ps2<0.0 && r>1.0){
					printf("screwy, psi^2=%g, r=%g, q=%g, ctheta=%g\n",ps2,r,q,ctheta);
					exit(1);
				}
				wfarray[iq][ir][ictheta]=ps2-1.0;
      }
    }
  }
}

double CKernelWF::GetPsiSquared(int iq,int ir,int ictheta){
  if(iq<nqmax && ir<nrmax && ictheta<=nctheta)
    return 1.0+wfarray[iq][ir][ictheta];
  else return 1.0;  
}

double CKernelWF::GetPsiSquared(int iq,int ir,double ctheta){
  double wa,wb;
  int ictheta;
  if(iq<nqmax && ir<nrmax){
    if(IDENTICAL) ctheta=fabs(ctheta);
    ictheta=lrint(floor((1.0-ctheta)/delctheta));
    wb=((1.0-ctheta)-ictheta*delctheta)/delctheta;
    wa=1.0-wb;
    return 1.0+wa*wfarray[iq][ir][ictheta]+wb*wfarray[iq][ir][ictheta+1];
  }
  else return 1.0;
}

double CKernelWF::GetPsiSquared(int iq,double r,double ctheta){
  int ir;
  double wa,wb;
  if(iq<nqmax){
    ir=lrint(floor((r-0.5*delr)/delr));
    if(ir<0) ir=0;
    wb=(r-(ir+0.5)*delr)/delr;
    wa=1.0-wb;
    return wa*GetPsiSquared(iq,ir,ctheta)+wb*GetPsiSquared(iq,ir+1,ctheta);
  }
  else return 1.0;
}

double CKernelWF::GetPsiSquared(double q,double r,double ctheta){
  int iq;
  double wa,wb;
  iq=lrint(floor((q-0.5*delq)/delq));
  if(iq<0) iq=0;
  wb=(q-(iq+0.5)*delq)/delq;
  wa=1.0-wb;
  return wa*GetPsiSquared(iq,r,ctheta)+wb*GetPsiSquared(iq+1,r,ctheta);
}

double CKernelWF::GetDELR(){
  return delr;
}

double CKernelWF::GetDELQ(){
  return delq;
}

double CKernelWF::GetDELCTHETA(){
  return delctheta;
}

int CKernelWF::GetNRMAX(){
  return nrmax;
}

int CKernelWF::GetNQMAX(){
  return nqmax;
}

int CKernelWF::GetNCTHETA(){
  return nctheta;
}

bool CKernelWF::GetIDENTICAL(){
  return IDENTICAL;
}

void CKernelWF::WriteData(string datadir){
  double q;//,r;
  int iq,ir,ictheta;
  double meanwf2;
  char filename[80];
  FILE *outfile;
  char shellcommand[120];
  //sprintf(shellcommand,"#!sh; if [ ! -e  %s ]; then mkdir -p %s; fi; \0",
  //  datadir,datadir);
  sprintf(shellcommand,"mkdir -p %s",datadir.c_str());
  system(shellcommand);
  for(iq=0;iq<nqmax;iq++){
    meanwf2=0.0;
    q=(0.5+iq)*delq;
    sprintf(filename,"%s/q%07.2f.tmp",datadir.c_str(),q);
    outfile=fopen(filename,"w");
    fprintf(outfile,"%d %d\n",nrmax,nctheta);
    for(ir=0;ir<nrmax;ir++){
      for(ictheta=0;ictheta<=nctheta;ictheta++){
				meanwf2+=wfarray[iq][ir][ictheta];
				fprintf(outfile,"%17.10e ",wfarray[iq][ir][ictheta]);
      }
      fprintf(outfile,"\n");
    }
    meanwf2=meanwf2/double(nctheta*nrmax);
    //printf("mean wfarray for iq=%d is %g\n",iq,meanwf2);
    fclose(outfile);
  } 
}

void CKernelWF::ReadData(string datadir){
  double q;
  double meanwf2;
  int iq,ir,ictheta,nrmaxread,ncthetaread;
  char filename[80];
  FILE *infile;
	
  for(iq=0;iq<nqmax;iq++){
    meanwf2=0.0;
    q=(0.5+iq)*delq;
    sprintf(filename,"%s/q%07.2f.tmp",datadir.c_str(),q);
    infile=fopen(filename,"r");
    fscanf(infile,"%d %d\n",&nrmaxread,&ncthetaread);
    if(nrmaxread!=nrmax || ncthetaread!=nctheta){
      printf("CKernelWF : trying to read file with wrong dimensions\n");
      printf("nrmaxread=%d, nrmax=%d, ncthetaread=%d, nctheta=%d\n",
				nrmaxread,nrmax,ncthetaread,nctheta);
      exit(1);
    }
    for(ir=0;ir<nrmax;ir++){
      for(ictheta=0;ictheta<=nctheta;ictheta++){
				fscanf(infile,"%lf ",&wfarray[iq][ir][ictheta]);
				meanwf2+=wfarray[iq][ir][ictheta];
      }
    }
    fclose(infile);
    meanwf2=meanwf2/double(nctheta*nrmax);
    //printf("mean wfarray for iq=%d is %g\n",iq,meanwf2);
  } 
}



