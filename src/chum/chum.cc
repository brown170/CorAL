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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <math.h>
#include "random.h"
#include "parametermap.h"
#include "cheezyparser.h"
#include "constants.h"
   
#define NTYPES 6

#undef USE_OLD_INPUT_SCHEME

using namespace std;

string flowprof("transverse");

void getCoreMom(double ans[4], const double &temp, const double &mass, const double vflow[4]);
void getCorePos(double ans[4], const double &Rx, const double &Ry, const double &Rz, const double &tau);
void getFlowVelocity(double ans[4], const double r[4], const double R, const double vmax=1./sqrt(3.));
void getHaloPos(double ans[4], const double r[4], const double p[4], const double &mass_w, const double &tau_w);
void getHaloMom(double ans[4], const double p_res[4], const double mass, const double mass_res);
bool passesCut( const double p[4], const double etamin, const double etamax, 
                const double ymin, const double ymax, const double pTmin, const double pTmax);
double dabs(const double& x){if (x>=0) {return x;} else {return -x;}}

// ------------- getHelp -------------
//! Print usage information & quit
void getHelp(void){
    cout<<"\nUsage: chum <options> [inputFile.dat]"<<endl;
    cout<<"    -h, -help, --help : Prints this message"<<endl;
    exit(0);
}

CRandom Randomizer(197);

/**-------------------------------------------------------
  *
  *  Core-Halo Ur-Model (CHUM)
  *  This is a simple implementation of the Core-Halo model
  *  based off of Scott-Pratt's PhaseMaker code, but with  
  *  several simple flow patterns implemented as well as 
  *  a rho-meson induced halo as well as an exponential 
  *  decaying core (which also gives a halo!)
  *
  *-------------------------------------------------------
  */
int main(int argc, char* argv[]){

    cout << "*** CHUM: the Core-Halo Ur-Model ***"<<endl;
    cout << endl;

    MESSAGE << CMessage::warning;

    // Parse command line
    if (argc==1) getHelp();
    string paramFile("");
    vector<string> modeList;
    for (int iarg = 1; iarg<argc; ++iarg){
        string sarg(argv[iarg]);
        if (sarg=="-help") getHelp();
        if (sarg=="--help") getHelp();
        if (sarg=="-h") getHelp();
        if (sarg.substr(0,1)=="-") modeList.push_back(sarg);
        else paramFile = sarg;
    }    
    
    // Initialize Variables
    double p[4],r[4],vflow[4];
    double mass[NTYPES]={938.3,939.6,139.58,139.58,134.9766,493.7};
    int ident[NTYPES]={2212,2112,211,-211,111,321};//p, n, pi+, pi-, K+/-
    double T, Rx, Ry, Rz, tau_fo, tau_w, fraction;
    double etamin,etamax,ymax,ymin,pTmax,pTmin;
    double R_Au, mass_w;
    int id,iblankline,iseed;
    string filename;
    int NUMPARTICLES = 100000;

    // Setting derived parameters
    tau_w = 23.0; // omega meson lifetime, fm/c
    R_Au  = 1.2*pow(197.,1./3.); // Au radius, fm
    id    = 2;   // pi+ only
    mass_w= 782.6; // mass of omega meson, MeV

    // setting fiducial cuts
    etamin=0.0;
    ymin=0.0;
    pTmin=0.0;
    etamax=200.0;
    ymax=200.0;
    pTmax=2e6;
    
    // Get user-setable model parameters
#ifdef USE_OLD_INPUT_SCHEME
    cout << "\n";
    cout << "What is the temperature? (in MeV)\n";
    cin >> T;
    cout << "What are Rx, Ry, Rz of the core and tau_fo of the rxn? (in fm)\n";
    cin >> Rx >> Ry >> Rz >> tau_fo;
    cout << "What fraction of particles come from the core? \n";
    cin >> fraction;
#else
    if (argc==1){
         cout<<"\nUsage: chum [input parameter file]"<<endl;
         exit(false);
    }
    ifstream pmapFile(argv[1]);
    if (!pmapFile.good()){
        cerr<< "Cannot open file "<<argv[1]<<endl;
        exit(-1);
    }
    parameterMap pmap;
    pmapFile>>pmap;
    T        = parameter::getD(pmap,"T",175.0);
    Rx       = parameter::getD(pmap,"Rx",4.0);
    Ry       = parameter::getD(pmap,"Ry",4.0);
    Rz       = parameter::getD(pmap,"Rz",4.0);
    tau_fo   = parameter::getD(pmap,"tau_fo",10.0);
    fraction = parameter::getD(pmap,"f",1.0);
    id       = parameter::getI(pmap,"id",id);
    flowprof = parameter::getS(pmap,"flowprofile",flowprof);
    etamin   = parameter::getD(pmap,"etamin",etamin);
    etamax   = parameter::getD(pmap,"etamax",etamax);
    ymin     = parameter::getD(pmap,"ymin",ymin);
    ymax     = parameter::getD(pmap,"ymax",ymax);
    pTmin    = parameter::getD(pmap,"pTmin",pTmin);
    pTmax    = parameter::getD(pmap,"pTmax",pTmax);
    filename = parameter::getS(pmap,"filename",(string)"thermal_gauss01.dat");
    iseed    = parameter::getI(pmap,"initial_seed",197);
    NUMPARTICLES = parameter::getI(pmap,"num_particles",100000);
#endif

    Randomizer.reset(iseed);
    
    cout << "Core radii are " << Rx << " fm, " << Ry << " fm, " << Rz << " fm\n";
    cout << "Core freezeout duration is " << tau_fo << " fm/c" <<endl;
    cout << "Core fraction is " << fraction << "\n";
    cout << "Particle type is " << ident[id] << " with mass " << mass[id] <<" MeVc^2"<<endl;
    cout << "Accepting particles w/ "<< etamin << " < eta < "<<etamax<<endl;
    cout << "Accepting particles w/ "<< ymin << " < y < "<<ymax<<endl;
    cout << "Accepting particles w/ "<< pTmin << " < pT < "<<pTmax<<" MeV"<<endl;
    cout << "Output filename "<< filename <<endl;
    cout << "Random seed "<< iseed <<endl;
    cout << "Number of particles generated "<< NUMPARTICLES <<endl;

    if      (flowprof==(string)"noflow")     cout << "Flow turned off"        << endl;
    else if (flowprof==(string)"radial")     cout << "Radial flow"            << endl;
    else if (flowprof==(string)"transverse") cout << "Transverse flow only"   << endl;
    else if (flowprof==(string)"bjorken")    cout << "Longitudinal flow only" << endl;
    else {
        cerr << "Illegal flow profile:" << flowprof << endl;
        exit(false);
    }
    // Main Loop over impact parameters
    ofstream fout(filename.c_str());
    for(iblankline=0;iblankline<3;iblankline++){
        fout << "line number " << iblankline << ", blah blah blah\n";
    }
    fout << 1 << " " 
         << NUMPARTICLES << " " 
         << setw(9) << setprecision(6) << 0.0 << " " 
         << setw(9) << setprecision(6) << 0.0 << "\n";
    // Loop over all particles
    int i=0;
    while (i<NUMPARTICLES) {
        getCorePos(r,Rx,Ry,Rz,tau_fo);
        getFlowVelocity(vflow,r,R_Au);
        if (Randomizer.ran()<=fraction){
            getCoreMom(p,T,mass[id],vflow);
        } else {
            getCoreMom(p,T,mass_w,vflow);
            getHaloPos(r,r,p,mass_w,tau_w);
            getHaloMom(p,p,mass[id],mass_w);
        }
        if (passesCut(p,etamin,etamax,ymin,ymax,pTmin,pTmax)){
            // 
            // output OSCAR formatted output
            //
            fout << i+1 << " "<< ident[id] << " "
                << setw(9) << p[1]/1000.0     << " "
                << setw(9) << p[2]/1000.0     << " "
                << setw(9) << p[3]/1000.0     << " "
                << setw(9) << p[0]/1000.0     << " "
                << setw(9) << mass[id]/1000.0 << " "
                << setw(9) << r[1]            << " "
                << setw(9) << r[2]            << " "
                << setw(9) << r[3]            << " "
                << setw(9) << r[0]            << "\n";
            i++;
        }
    }
    return 0;
}

//**********************************************************
//
//  Generates a momentum from the core, sampling  
//  p from exp(-p.vflow/T+m/T)
//
//**********************************************************
void getCoreMom(double ans[4], const double &temp, const double &mass, const double vflow[4])
{
    double Tx(0.), x(0.), mu(0.), gammainv2(1.), Txm(0.);
    double vmag = sqrt(vflow[1]*vflow[1]+vflow[2]*vflow[2]+vflow[3]*vflow[3]);
    x = Randomizer.ran_exp();
    mu=2.*Randomizer.ran()-1.; // will be angle between p and vflow
    Tx = temp*x;
    Txm=Tx+mass;
    gammainv2 = 1.0-vmag*vmag*mu*mu;
    double pmag=(vmag*mu*Txm+sqrt(Txm*Txm-gammainv2*mass*mass))/gammainv2;
    // get rest of spherical coords
    double sthet=sqrt(1.0-mu*mu);
    double phi=2.0*PI*Randomizer.ran();
    // compute 4-mom
    ans[0]=sqrt(mass*mass+pmag*pmag);
    ans[1]=pmag*sthet*cos(phi);
    ans[2]=pmag*sthet*sin(phi);
    ans[3]=pmag*mu;                          
}
//**********************************************************
//
//  Generates a position in the core
//
//**********************************************************
void getCorePos(double ans[4], const double &Rx, const double &Ry, const double &Rz, const double &tau)
{
    double tmp;
    ans[0] = Randomizer.ran_exp();
    Randomizer.gauss2(&(ans[1]),&tmp);
    Randomizer.gauss2(&(ans[2]),&(ans[3]));
    ans[0]=tau*ans[0];
    ans[1]=Rx*ans[1];
    ans[2]=Ry*ans[2];
    ans[3]=Rz*ans[3];
}
//**********************************************************
//
//  compute flow velocity 
//
//**********************************************************
void getFlowVelocity(double ans[4], const double r[4], const double R, const double vmax)
{
    ans[0]=1.0; for (int i=1;i<4;++i)ans[i]=0.0;
    double alpha = vmax/R; // flow param, c/fm

    if (flowprof=="radial") {
//      assume radial flow within Zajc's radial flow
//      model, but with maximal vflow of vmax
        double rmag = sqrt(r[1]*r[1]+r[2]*r[2]+r[3]*r[3]);
        if (rmag>=R) { for (int i=1;i<4;++i){ans[i] = vmax;} }
        else{ for (int i=1;i<4;++i){ans[i] = alpha*r[i];} }
    }
    if (flowprof=="transverse") {
//      assume transverse flow within simplified Zajc radial
//      flow model, but with maximal vflow of vmax
        double rT = sqrt(r[1]*r[1]+r[2]*r[2]);
        if (rT>=R) { for (int i=1;i<3;++i){ans[i] = vmax;} }
        else{ for (int i=1;i<3;++i){ans[i] = alpha*r[i];} }
        ans[3]=0.0; // no logitudinal flow
    }
    if (flowprof=="bjorken") {
//      cheezy version of Bjorken flow, essentially Zajc
//      model in Z (beam) direction only.  maximal flow of vmax
        if (r[3]>=R) { ans[3] = vmax;}
        else { ans[3] = alpha*r[3]; }
        ans[1]=0.0; ans[2]=0.0; // no transverse flow
    }

}
//**********************************************************
//
//  Generates a position in the halo
//
//**********************************************************
void getHaloPos(double ans[4], const double r[4], const double p[4], const double &mass_w, const double &tau_w)
{
    double t = Randomizer.ran_exp();
    for (int i=0;i<4;++i){ans[i]=r[i]+p[i]*tau_w*t/p[0];}
}
//**********************************************************
//
//  Generates a momentum in the halo
//
//**********************************************************
void getHaloMom(double ans[4], const double p_res[4], const double mass, const double mass_res)
{
    // Assume's 3-body decay, but pi+, pi-, and pi0 all have same mass
    int iparticle(static_cast<int>(10*Randomizer.ran()/3));
    
//    double E = p_res[0];
    double M = mass_res;
    double M1 = mass;
    double M2 = mass;
    double M3 = mass;
    double ES1, ES2, P1, P2, Cos12, Z;
      
    do {
	    // Generate E1 and E2 with the Momnte-Carlo method
	    do {
	        ES1 = Randomizer.ran() * (M - M2 - M3 - M1) + M1;
	        ES2 = Randomizer.ran() * (M - M1 - M3 - M2) + M2;
	    } while (ES1+ES2 > M); // The sum of both energies must be smaller than the resonance mass
      
	    P1  = sqrt(ES1*ES1 - M1*M1);
	    P2  = sqrt(ES2*ES2 - M2*M2);
	    Z = M - ES1 - ES2;
	    Z *= Z;
	    Cos12 = (Z - P1*P1 - P2*P2 - M3*M3)/(2*P1*P2);
	} while ((Cos12 < -1.0) || (Cos12 > 1.0)); // Cos Theta must exist (be within -1.0 to 1.0 )
          
    double Pxr2 = P2 * sqrt(1-Cos12*Cos12);
    double Pzr2 = P2*Cos12;
    double Pxr3 = - Pxr2;
    double Pzr3 = - (P1 + Pzr2);
    double P3 = sqrt(Pxr3*Pxr3+Pzr3*Pzr3);
    double ES3 = sqrt(M3*M3+P3*P3);

    // Generating Euler angles
    double Phi = Randomizer.ran()   * 2 * PI;
    double Ksi = Randomizer.ran()   * 2 * PI;
    double CosTh = Randomizer.ran() * 2.0 - 1.0;

    double sp = sin(Phi);
    double cp = cos(Phi);
    double sk = sin(Ksi);
    double ck = cos(Ksi);
    double st = sqrt(1.0-CosTh*CosTh);
    double ct = CosTh;

    // Rotating the whole system
    double Pxp1 = - st*ck * P1;
    double Pyp1 = st*sk * P1;
    double Pzp1 = ct * P1;
    
    double Pxp2 = (cp*ct*ck - sp*sk)  * Pxr2 + (-st*ck) * Pzr2;
    double Pyp2 = (-cp*ct*sk - sp*ck) * Pxr2 + (st*sk)  * Pzr2;
    double Pzp2 = cp*st               * Pxr2 + ct       * Pzr2;

    double Pxp3 = (cp*ct*ck - sp*sk)  * Pxr3 + (-st*ck) * Pzr3;
    double Pyp3 = (-cp*ct*sk - sp*ck) * Pxr3 + (st*sk)  * Pzr3;
    double Pzp3 = cp*st               * Pxr3 + ct       * Pzr3;
    
    double Vx = p_res[1]/p_res[0]; 
    double Vy = p_res[2]/p_res[0]; 
    double Vz = p_res[3]/p_res[0]; 
    
    ES1 = sqrt(M1*M1+Pxp1*Pxp1+Pyp1*Pyp1+Pzp1*Pzp1);
    ES2 = sqrt(M2*M2+Pxp2*Pxp2+Pyp2*Pyp2+Pzp2*Pzp2);
    ES3 = sqrt(M3*M3+Pxp3*Pxp3+Pyp3*Pyp3+Pzp3*Pzp3);

    double V2 = Vx*Vx + Vy*Vy + Vz*Vz;
    double Gamma = 1.0/sqrt(1-V2);
    
    // Boosting by the parent velocity
    double VP = Vx*Pxp1 + Vy*Pyp1 + Vz*Pzp1;
    double gvp = (Gamma - 1.0) * (1.0/V2) * VP;

    double Px1 = Pxp1 + (gvp + Gamma * ES1) * Vx;
    double Py1 = Pyp1 + (gvp + Gamma * ES1) * Vy;
    double Pz1 = Pzp1 + (gvp + Gamma * ES1) * Vz;
  
    VP = Vx*Pxp2 + Vy*Pyp2 + Vz*Pzp2;
    gvp = (Gamma - 1.0) * (1.0/V2) * VP;

    double Px2 = Pxp2 + (gvp + Gamma * ES2) * Vx;
    double Py2 = Pyp2 + (gvp + Gamma * ES2) * Vy;
    double Pz2 = Pzp2 + (gvp + Gamma * ES2) * Vz;
  
    VP = Vx*Pxp3 + Vy*Pyp3 + Vz*Pzp3;
    gvp = (Gamma - 1.0) * (1.0/V2) * VP;

    double Px3 =   Pxp3 + (gvp + Gamma * ES3) * Vx;
    double Py3 =   Pyp3 + (gvp + Gamma * ES3) * Vy;
    double Pz3 =   Pzp3 + (gvp + Gamma * ES3) * Vz;
    
    ES1 = sqrt(M1*M1+Px1*Px1+Py1*Py1+Pz1*Pz1);
    ES2 = sqrt(M2*M2+Px2*Px2+Py2*Py2+Pz2*Pz2);
    ES3 = sqrt(M3*M3+Px3*Px3+Py3*Py3+Pz3*Pz3);

    if (iparticle==1) {
        ans[0]=ES1;ans[1]=Px1;ans[2]=Py1;ans[3]=Pz1;
    } else if (iparticle==2) {
        ans[0]=ES2;ans[1]=Px2;ans[2]=Py2;ans[3]=Pz2;
    } else {
        ans[0]=ES3;ans[1]=Px3;ans[2]=Py3;ans[3]=Pz3;
    } 

}
//**********************************************************
//
//  Checks if particle momentum passes acceptance cuts
//
//**********************************************************
bool passesCut( const double p[4], const double etamin, const double etamax, const double ymin, const double ymax, const double pTmin, const double pTmax)
{
    double y=0.5*log((p[0]+p[3])/(p[0]-p[3]));
    double pT=sqrt(p[1]*p[1]+p[2]*p[2]);
// don't do psuedorapidity until have numerically stable implementation
//    double theta=atan(pT/p[3]);
//    double eta=-log(abs(tan(theta/2.0)));
    bool goody  = (dabs(y)>=ymin)&&(dabs(y)<=ymax);
    bool goodpT = (pT>=pTmin)&&(pT<=pTmax);
    bool goodeta= (dabs(y)>=etamin)&&(dabs(y)<=etamax);
//    if (goody&&goodpT&&goodeta)
//    cout << "y:"<< y <<" "<< ymin << " "<<ymax <<" "<< goody << " eta:"<<etamin << " "<<etamax  << " " << goodeta << 
//          " pT:"<<pT<<" "<<pTmin<<" " << pTmax << " " << goodpT<<" p:"<<p[0] << " "<<p[1]<< " "<<p[2]<<" "<< p[3]<<" "<<endl;
    return goody&&goodpT&&goodeta;
}

