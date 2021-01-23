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
#ifndef __INCLUDE_SOURCECALC_H
#define __INCLUDE_SOURCECALC_H
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include "arrays.h"
#include "parametermap.h"
#include "random.h"
#include "constants.h"

using namespace std;

class CSourceCalc{
 public:
  parameterMap spars;
  virtual void CalcS(CCHArray *A);
  virtual void CalcS(int lx, int ly, int lz, CCHArray *A);
  virtual void CalcS(CMCList *&lista, CMCList *&listb);
	virtual void CalcS(C3DArray *threed);
	virtual void GaussCFCalc(C3DArray *cf3d);

  void CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A);
  void CombineMCLists(CMCList *lista,CMCList *listb,CCHArray *A,int NMC);
  void CombineMCLists(CMCList *lista,CMCList *listb,C3DArray *threed);
  void ReadSPars(char *sparsfilename);
  void NormCheck(CCHArray *A);
  void NormCheck(C3DArray *threed);
  void CalcEffGaussPars(CCHArray *A);
  void CalcEffGaussPars(CCHArray *A,double &Rx,double &Ry,double &Rz,
			double &Xoff,double &Yoff,double &Zoff);
  CSourceCalc();
	CSourceCalc(string sparsfilename);
	CRandom *randy;
	virtual ~CSourceCalc(){};
private:
	
};

class CSourceCalc_GX1D : public CSourceCalc{
/*
  S = { N_G * lambda_G * exp(-r^2/4R^2)
        + N_X * lambda_X * exp(-[r^2/X^2+4R^4/X^4]) } * [r^2/(r^2+a^2)]^(L/2)
	where N_x and N_g are constants that make individual contributions
        integrate to unity for L=0
*/
 public:
  CSourceCalc_GX1D();
  void InitSPars();
  void SetSPars(double lambda,double Xfrac,double R,double X,double a);
  void CalcS(int lx,int ly,int lz,CCHArray *A);
  void CalcS(CCHArray *A);
};

class CSourceCalc_Gaussian : public CSourceCalc{
 public:
  CSourceCalc_Gaussian();
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset);
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset);
  void SetSPars(double lambdaset,double Rxset,double Ryset,double Rzset,
		double Xoffset,double Yoffset,double Zoffset,
		double Euler_phiset,double Euler_thetaset,
		double Euler_psiset);
  void CalcS(CCHArray *A);
  void CalcS(C3DArray *threed);
	void GaussCFCalc(C3DArray *cf3d);
  
 private:
  void CalcAlpha(double **alpha,CCHArray *A);
  void InitSPars();
  // S ~ exp{-alpha_ij (x_i-off_i)(x_j-off_j)} 
  // Calcs ignore zeroth components
};

class CSourceCalc_EllipticBlast : public CSourceCalc{
 public:
  CSourceCalc_EllipticBlast();
  void CalcS(CCHArray *A);
  void SetSPars(double Rxset,double Ryset,
		double Tauset,
		double BetaXset,double BetaYset,
		double Tset,double Ptset,
		double Phiset,double EtaGset,
		double Maset,double Mbset);
  void SetSPars(double Rset,double Tauset,
		double Betaset,double Tset,double Ptset);
 private:
  void Get_r(double *p,int nsample,double **r);
  void InitSPars();
};

class CSourceCalc_Blast : public CSourceCalc{
 public:
  CSourceCalc_Blast();
  void CalcS(CCHArray *A);
  void CalcS(CMCList *lista,CMCList *listb);
  void SetSPars(double lambdaset,
		double Rset,double Tauset,double DelTauset,
		double Betaset,double Tset,double Ptset,
		double EtaGset,double Maset,double Mbset);
  void SetSPars(double lambdaset,double Rset,double Tauset,double DelTauset);
 private:
  void GetMCList(double *p,CMCList *mclist);
  void InitSPars();
  double GetTau(double tau0,double deltau);
};

class CSourceCalc_OSCAR : public CSourceCalc{
 public:
  CSourceCalc_OSCAR();
  CSourceCalc_OSCAR(string sparsfilename);
  //void CalcS(CCHArray *A);
	void CalcS(CMCList *&lista,CMCList *&listb);
	void SetSPars(double Pt_set,double delpt_set,double phimin_set,double phimax_set,double ymin_set,double ymax_set);
	void SetIDs(int *idlista,int nida,int *idlistb,int nidb);
private:
	int *idlist_a,*idlist_b,nid_a,nid_b;
  void ReadR(double **ra,int &na,double **rb,int &nb);
  void InitSPars();
	bool IDMatch(int ident,int *idlist,int nid);
	bool Check(double *p,double *r,double m,double **ra,int &n);
};

#endif
