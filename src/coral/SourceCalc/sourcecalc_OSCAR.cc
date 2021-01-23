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
#ifndef __INCLUDE_SOURCECALC_OSCAR
#define __INCLUDE_SOURCECALC_OSCAR
#include "sourcecalc.h"

using namespace std;

CSourceCalc_OSCAR::CSourceCalc_OSCAR(){
	InitSPars();	
	parameter::PrintPars(spars);
	randy=new CRandom(1234);
}

CSourceCalc_OSCAR::CSourceCalc_OSCAR(string sparsfilename){
	InitSPars();
	parameter::ReadParsFromFile(spars,sparsfilename);
	parameter::PrintPars(spars);
	randy=new CRandom(1234);
}

void CSourceCalc_OSCAR::InitSPars(){
	// DEFAULT VALUES
	parameter::set(spars,"PT",600.0);
	parameter::set(spars,"DELPT",20.0);
	parameter::set(spars,"PHIMIN_DEG",0.0);
	parameter::set(spars,"PHIMAX_DEG",360.0);
	parameter::set(spars,"YMIN",-1.0);
	parameter::set(spars,"YMAX",1.0);
	parameter::set(spars,"NMAX",20000);
	parameter::set(spars,"OSCARfilename",string("UNDEFINED"));
	parameter::set(spars,"NEVENTSMAX",10000);
	parameter::set(spars,"ETA_GAUSS",1.2);
}

void CSourceCalc_OSCAR::SetSPars(double PT_set,double DELPT_set,double PHIMIN_DEG_set,double PHIMAX_DEG_set,double YMIN_set,double YMAX_set){
	parameter::set(spars,"PT",PT_set);
	parameter::set(spars,"DELPT",DELPT_set);
	parameter::set(spars,"PHIMIN_DEG",PHIMIN_DEG_set);
	parameter::set(spars,"PHIMAX_DEG",PHIMAX_DEG_set);
	parameter::set(spars,"YMIN",YMIN_set);
	parameter::set(spars,"YMAX",YMAX_set);
}

void CSourceCalc_OSCAR::SetIDs(int *idlista,int nida,int *idlistb,int nidb){
	idlist_a=idlista;
	nid_a=nida;
	nid_b=nidb;
	idlist_b=idlistb;
}

void CSourceCalc_OSCAR::CalcS(CMCList *&lista, CMCList *&listb){
	double **ra,**rb;
	int ia,ib,na=0,nb=0;
	bool AEQUALB=parameter::getB(spars,"AEQUALB",true);
	int NMAX=parameter::getI(spars,"NMAX",50000);

	ra=new double *[NMAX];
	for(ia=0;ia<NMAX;ia++) ra[ia]=new double[4];
	if(AEQUALB){
		rb=ra;
	}
	else{
		rb=new double *[NMAX];
		for(ib=0;ib<NMAX;ib++) rb[ib]=new double[4];
	}
	ReadR(ra,na,rb,nb);

	if(lista==NULL) lista=new CMCList(na);
	else if(lista->GetNMC()!=na) lista->Resize(na);
	if(!AEQUALB){
		if(listb==NULL) listb=new CMCList(nb);
		else if(listb->GetNMC()!=nb) listb->Resize(nb);
	}

	for(int i=0;i<na;i++) lista->SetR(i,ra[i]);
	if(lista!=listb){
		for(int i=0;i<nb;i++) listb->SetR(i,rb[i]);
	}

	for(ia=0;ia<NMAX;ia++){
		delete [] ra[ia];
	}
	delete [] ra;
	if(!AEQUALB){
		for(ib=0;ib<NMAX;ib++){
			delete [] rb[ib];
		}
		delete [] rb;
	}
	printf("______ FINISHED CREATING MCLISTS ___________\n");

}

void CSourceCalc_OSCAR::ReadR(double **ra, int &na, double **rb, int &nb){
	FILE *oscarfile;
	double r[4],p[4];
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double mass,rdummy1,rdummy2;
	int ident,idummy,i,alpha;
	int ndummy_header=3,ndummy_betweenevents=0;
	int npart,npartmax,ipart,ievent;
	int NEVENTSMAX=parameter::getI(spars,"NEVENTSMAX",100);
	string OSCARfilename=parameter::getS(spars,"OSCARfilename","UNDEFINED");
	bool AEQUALB=parameter::getB(spars,"AEQUALB",false);
	char dummy[160];
	long long int nsigmay=0;
	//double y,eta,sigmay=0.0,sigmaz=0.0,z;
	//double taucompare=parameter::getD(spars,"TAUCOMPARE",12.0);
	//taucompare=12.0;
	na=nb=0;
	printf("Opening %s\n",OSCARfilename.c_str());
	oscarfile=fopen(OSCARfilename.c_str(),"r");
	// Read Header Info
	for(i=0;i<ndummy_header;i++){
		fgets(dummy,120,oscarfile);
		//printf("dummy line=%s\n",dummy);
	}

	//for(ievent=1;ievent<=nevents;ievent++){
	ievent=0;
	do{
		ievent+=1;
		fscanf(oscarfile,"%d %d %lf %lf",&idummy,&npartmax,&rdummy1,&rdummy2);
		fgets(dummy,120,oscarfile);
		//printf("ievent=%d, idummy=%d, npartmax=%d\n",ievent,idummy,npartmax);
		for(npart=0;npart<npartmax;npart++){
			fscanf(oscarfile,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf",
				&ipart,&ident,&p[1],&p[2],&p[3],&p[0],&mass,&r[1],&r[2],&r[3],&r[0]);
			for(alpha=0;alpha<4;alpha++) p[alpha]*=1000.0;
			/* SET KINEMATIC VARIABLES */

			/*
			if(r[0]<40){
				eta=atanh(r[3]/r[0]);
				y=atanh(p[3]/p[0]);
				nsigmay+=1;
				sigmay+=(eta-y)*(eta-y);
				z=sqrt(r[0]*r[0]-r[3]*r[3])*sinh(eta-y);
				sigmaz+=z*z;
			}
			*/
			if(IDMatch(ident,idlist_a,nid_a)) Check(p,r,MA,ra,na);
			if(!AEQUALB) if(IDMatch(ident,idlist_b,nid_b)) Check(p,r,MB,rb,nb);
		}
		for(i=0;i<ndummy_betweenevents;i++){
			if(!feof(oscarfile)) fgets(dummy,120,oscarfile);
			//printf("dummy line=%s\n",dummy);
		}
	} while(ievent<NEVENTSMAX && !feof(oscarfile));
	fclose(oscarfile);
	if(AEQUALB) nb=na;
	//sigmay=sqrt(sigmay/nsigmay);
	//sigmaz=sqrt(sigmaz/nsigmay);
	//printf("OSCAR file read: %d events, na=%d, nb=%d, sigma_y=%g, sigma_z=%g\n",ievent,na,nb,sigmay,sigmaz);
	printf("OSCAR file read: %d events, na=%d, nb=%d\n",ievent,na,nb);
}

bool CSourceCalc_OSCAR::Check(double *p, double *r, double m, double **ra, int &n){
	double PT=parameter::getD(spars,"PT",600.0);
	double DELPT=parameter::getD(spars,"DELPT",20.0);
	double YMIN=parameter::getD(spars,"YMIN",-1.0);
	double YMAX=parameter::getD(spars,"YMAX",1.0);
	double PHIMIN=2.0*PI*parameter::getD(spars,"PHIMIN_DEG",0.0)/360.0;
	double PHIMAX=2.0*PI*parameter::getD(spars,"PHIMAX_DEG",0.0)/360.0;
	double ETA_GAUSS=parameter::getD(spars,"ETA_GAUSS",1.2);
	double MA=parameter::getD(spars,"MA",139.57);
	double MB=parameter::getD(spars,"MB",139.57);
	double pt,phi,eta;
	double gammav=PT/(MA+MB);
	double gamma=sqrt(1.0+gammav*gammav);
	double pttarget=PT*m/(MA+MB);
	double rout,rlong,rside,sinhy,coshy,tau,vperp,y;
	const double TAUCOMPARE=12.0;
	int NMAX=parameter::getI(spars,"NMAX",50000);
	bool success=false;
	
	pt=sqrt(p[1]*p[1]+p[2]*p[2]);
	if(fabs(pt-pttarget)<DELPT){
		if(p[1]!=p[1] || p[2]!=p[2] || p[3]!=p[3]){
			printf("bad particle has nan, p=(%g,%g,%g)\n",p[1],p[2],p[3]);
			return false;
		}
		phi=atan2(p[1],p[2]);
		if( ((PHIMIN<PHIMAX)&&(phi>PHIMIN && phi<PHIMAX) )
		|| ((PHIMIN>PHIMAX)&&(phi>PHIMIN||phi<PHIMAX)) ){
			y=atanh(p[3]/p[0]);
			if(y>YMIN && y<YMAX){
				p[0]=sqrt(pt*pt+p[3]*p[3]+m*m);
				rout=(p[1]*r[1]+p[2]*r[2])/pt;
				rside=(p[1]*r[2]-p[2]*r[1])/pt;
				sinhy=sinh(y);
				coshy=cosh(y);
				rlong=coshy*r[3]-sinhy*r[0];
				tau=coshy*r[0]-sinhy*r[3];
				eta=asinh(rlong/tau);
				if(randy->ran()<exp(-0.5*eta*eta/(ETA_GAUSS*ETA_GAUSS))){
				//printf("%g %g %g %g\n",tau,rout,rside,rlong);
					vperp=pt/sqrt(m*m+pt*pt);
					rout=rout-vperp*(tau-TAUCOMPARE);
					tau=TAUCOMPARE;
					rout=gamma*rout-gammav*tau;
					ra[n][0]=0.0;
					ra[n][1]=rout;
					ra[n][2]=rside;
					ra[n][3]=rlong;
					n+=1;
					success=true;
					if(n==NMAX){
						printf("TOO MANY PARTICLES FIT CRITERIA, increase parameter NMAX=%d\n",NMAX);
						exit(1);
					}
				}
			}
		}
	}
	return success;
}

bool CSourceCalc_OSCAR::IDMatch(int ident, int *idlist, int nid){
	int i;
	bool answer=false;
	for(i=0;i<nid;i++) if(ident==idlist[i]) answer=true;
	return answer;
}

#endif
