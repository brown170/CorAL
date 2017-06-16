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
#include "oscar.h"
#include "constants.h"
#include "misc.h"
#include "message.h"
#include <fstream>

// ------------------ readOSCARFile ------------------
//! reads in OSCAR formatted file of last interactions phase-space points and stores them in a vector
vector<COSCARLine> readOSCARFile(string oscarFile){
    if ( !file_exists( oscarFile ) ) 
        throw MESSAGE << "OSCAR file "<< oscarFile <<" not found. Halting!!"<<ENDM_FATAL;
    vector<COSCARLine> result;
    ifstream iFile(oscarFile.c_str());
    char dummy_header[80],endline_tester;
    cout << "Reading in "<<oscarFile<<"\n";
    cout << "    " << "Header info:\n";
    for(int i=0;i<3;++i)
    {
        iFile.get(dummy_header,80,'\n');
        if (iFile.get(endline_tester) && endline_tester!='\n') 
            MESSAGE<< "Line a little too long, so read-in screwed up.  Halting. \n"<<ENDM_FATAL;
        cout << "        " << dummy_header<<"\n";
    }  
    int num_particles,n1;
    {
        double dummy1, dummy2;
        // set the number of particles in the freeze-out distribution  
        iFile >> n1 >> num_particles >> dummy1 >> dummy2;
        cout << "       Event " << n1 << " has " << num_particles << " particles. " << dummy1 << " " << dummy2 << "\n";
    }
    long lineCount = 0;
    long totalLines=0;
    while (iFile.good()) {
        lineCount++;
	if(lineCount > num_particles){
	  double dummy1, dummy2;
	  // set the number of particles in the freeze-out distribution  
	  iFile >> n1 >> num_particles >> dummy1 >> dummy2;
	  cout << "       Event " << n1 << " has " << num_particles << " particles. " << dummy1 << " " << dummy2 << "\n";
	  lineCount = 0;
	} else {
	  COSCARLine newLine;
	  iFile>>newLine;
	  result.push_back(newLine);
	  totalLines++;
	}
    }
    cout<<"OSCAR File had "<<totalLines<<" particles in "<<n1<<" events."<<endl;
    return result;
}

// ------------------ filterPID ------------------
//! keep only particles from a vector of particles which have particle ID's that match that requested
vector<COSCARLine> filterPID(const vector<COSCARLine> &lines, int pid){
  vector<COSCARLine> tLines;
  for (vector<COSCARLine>::const_iterator itr=lines.begin(); itr!=lines.end(); itr++){
    if (pid==itr->pid) tLines.push_back(*itr);
  }
  return tLines;
}

// ------------------ getOSCARSource1d ------------------ 
//! takes the particles from a vector of OSCAR-type particles and forms the 1d (angle-averaged) two-particle imaged source function
CSourceFtn1dHisto getOSCARSource1d(vector<COSCARLine> lines, parameterMap p){

    // Initialization
    CSourceFtn1dHisto result;
    parameterMap souMap; 
    if (p.find("source_settings")!=p.end()) souMap = parameter::getMap( p,"source_settings" );
    else souMap = p;
    result.Read(souMap);
    double q_cut_oscar   = parameter::getD(p,"q_cut_oscar",0.060);
    int max_number_pairs = parameter::getI(p,"max_number_pairs",1000000);
    int totalPairs       = 0; // the full sum of the unnorm'ed source

    // Define tmp variables 
    double q_lab[4],q_cm[4],P[4],beta[4],sep_lab[4],sep_cm[4];
    double sqr_beta(0.), sqr_q_cm(0.);
    bool good_boost;
    double rsep_cm;

    cout<<"    Zeroing 1D source array\n";
    for (int i=0;i<result.ndata;++i) {result.data[i]=0;}

    cout<<"    Constructing 1D source\n";
    // loop over pairs
    for (vector<COSCARLine>::iterator p1=lines.begin(); p1!=lines.end(); ++p1){
        vector<COSCARLine>::iterator p2=p1;
        ++p2;
        while ( (totalPairs<max_number_pairs) && (p2!=lines.end()) ) {
            for (int i=0;i<4;++i){
                q_lab[i]=0.5*(p1->p[i]-p2->p[i]);  // relative mom. is 1/2 the mom. diff. of pair
                P[i]=p1->p[i]+p2->p[i];            // total mom. is total of pair mom.
                beta[i]=P[i]/P[0];
                sep_lab[i]=(p1->x[i]-p2->x[i]);
            }
            Misc::lorentz(beta,q_lab,q_cm);
            sqr_beta=beta[0]*beta[0];
            sqr_q_cm=q_cm[0]*q_cm[0];
            for (int i=1; i<4; ++i) {
                sqr_beta -= beta[i]*beta[i];
                sqr_q_cm -= q_cm[i]*q_cm[i];
            }
            good_boost=(sqr_beta >= 0.);
            // Check for low relative momentum q<q_cut 
            if ((-sqr_q_cm<q_cut_oscar) && good_boost ){  
                totalPairs+=1;
                // This is a good pair, let's bin it up the source function(s)...
                Misc::lorentz(beta,sep_lab,sep_cm); //get pair separation in CM frame
                // bin up 1D source
                rsep_cm=sqrt(sep_cm[1]*sep_cm[1]+sep_cm[2]*sep_cm[2]+sep_cm[3]*sep_cm[3]);
                int iBin = result.whatBin(rsep_cm,false);
                if (iBin>=0 && iBin<result.ndata) result.data[iBin]+=1;
            }
            ++p2;
        }
    }

    cout<<"    Normalizing 1D source\n";
    double numPairs1D,rweight,rval,dr;
    for (int i=0;i<result.ndata;++i) {
        rval    = result.midBin(i);
        dr      = result.binWidth(i);
        rweight = (rval*rval+dr*dr/12.0)*dr*4.0*PI;
        numPairs1D = result.data[i]; 
        result.data[i]=result.data[i]/rweight/totalPairs;
        if (result.data[i]>0.) result.uncert[i]=result.data[i]
            *sqrt(1.0/double(numPairs1D)+1.0/double(numPairs1D)/double(numPairs1D)/double(totalPairs));
    }

    cout << "    Cross-checking angle averaged source function integral...\n";
    double int_sou=0.0,int_sou2=0.0;
    for(int i=0;i<result.ndata;++i){
        int_sou  += 4.0*PI*result.midBin(i)*result.midBin(i)*result.binWidth(i)*result.data[i];
        int_sou2 += (result.midBin(i)*result.midBin(i)+result.binWidth(i)*result.binWidth(i)/12.0)*result.binWidth(i)*4.0*PI*result.data[i];
    }
    cout << "        Naive integral of source: " << int_sou << "\n";
    cout << "        Not so naive integral of source: " << int_sou2 << "\n";
    cout << "        Number of good pairs = "<<totalPairs<<"\n";
    return result;
}

// ------------------ getOSCARSource3d ------------------ 
//! takes the particles from a vector of OSCAR-type particles and forms the 3d two-particle imaged source function
CSourceFtn3dHisto getOSCARSource3d(vector<COSCARLine> lines, parameterMap p){

    // Initialization
    CSourceFtn3dHisto result;
    result.Read(p);
    double q_cut_oscar   = parameter::getD(p,"q_cut_oscar",0.060);
    int max_number_pairs = parameter::getI(p,"max_number_pairs",1000000);
    int totalPairs       = 0; // the full sum of the unnorm'ed source

    // Determine which of Side, Out and Long is x, y, and z
    int iSide = parameter::getI(p,"cartesian_side_index",0);
    int iOut  = parameter::getI(p,"cartesian_out_index",1);
    int iLong = parameter::getI(p,"cartesian_long_index",2);
    
    // Define tmp variables 
    double r[3], S[3], O[3], L[3];
    L[0]=0.0;L[1]=0.0;L[2]=1.0;
    double q_lab[4],q_cm[4],P[4],beta[4],sep_lab[4],sep_cm[4];
    double sqr_beta(0.), sqr_q_cm(0.);
    bool good_boost;
    double rS,rO,rL,maxS,minS,maxO,minO,maxL,minL;

    cout<<"Zeroing 3D source array\n";
    for (int i=0;i<result.ndata;++i) {result.data[i]=0;}

    cout<<"Constructing 3D source\n";
    // loop over pairs
    for (vector<COSCARLine>::iterator p1=lines.begin(); p1!=lines.end(); ++p1){
        vector<COSCARLine>::iterator p2=p1;
        while ( (totalPairs<max_number_pairs) && (p2!=lines.end()) ) {
            p2++;
            for (int i=0;i<4;++i){
                q_lab[i]=0.5*(p1->p[i]-p2->p[i]);  // relative mom. is 1/2 the mom. diff. of pair
                P[i]=p1->p[i]+p2->p[i];            // total mom. is total of pair mom.
                beta[i]=P[i]/P[0];
                sep_lab[i]=(p1->x[i]-p2->x[i]);
            }
            Misc::lorentz(beta,q_lab,q_cm);
            sqr_beta=beta[0]*beta[0];
            sqr_q_cm=q_cm[0]*q_cm[0];
            for (int i=1; i<4; ++i) {
                sqr_beta -= beta[i]*beta[i];
                sqr_q_cm -= q_cm[i]*q_cm[i];
            }
            good_boost=(sqr_beta >= 0.);
            // Check for low relative momentum q<q_cut 
            if ((-sqr_q_cm<q_cut_oscar) && good_boost ){  
                totalPairs+=1;
                // This is a good pair, let's bin it up the source function(s)...
                Misc::lorentz(beta,sep_lab,sep_cm); //get pair separation in CM frame
                // bin up 3D source
                // define unit three-vectors of Bertsch-Pratt coords.
                // O is unit vector || to the part of P that is perp to L 
                O[0]=P[1];O[1]=P[2];O[2]=0.0;
                double Omag = sqrt(O[0]*O[0]+O[1]*O[1]);
                O[0]/=Omag;O[1]/=Omag;
                // S is perp to both O and L
                S[0]=O[1];S[1]=-O[0];S[2]=0.0;           
                // break sep_cm into SOL coords and figure out 
                // the bins they correspond too
                rS=0.0;rO=0.0;rL=0.0;
                for (int i=0;i<3;++i) {
                    rS+=sep_cm[i+1]*S[i];
                    rO+=sep_cm[i+1]*O[i];
                    rL+=sep_cm[i+1]*L[i];
                }
                // Set x, y, z according to specified coordinate system
                r[iSide] = rS;
                r[iOut]  = rO;
                r[iLong] = rL;
                // Generate the unnormalized source histogram
                result.data[result.whatBin(r[0],r[1],r[2])]+=1;
                // Get some stats on the r_cm's
                maxS=std::max(rS,maxS);
                minS=std::min(rS,minS);
                maxO=std::max(rO,maxO);
                minO=std::min(rO,minO);
                maxL=std::max(rL,maxL);
                minL=std::min(rL,minL);
            }
        }
        break;
    }

    cout << "Normalizing 3D source function...\n";
    double rweight=result.binVolume();
    double npairs;
    for (int i=0;i<result.ndata;++i) {
        result.uncert[i]=0.0;
        if (result.data[i]>0) {
            npairs=result.data[i];
            result.data[i]/=double(totalPairs)*rweight;
            result.uncert[i]=result.data[i]*sqrt(1.0/npairs+1.0/npairs/npairs/double(totalPairs));
        }
    }
    cout << "Some parameters from the legal pairs:\n";
    cout << "  maxS:" << maxS << ", minS:" << minS << "\n";
    cout << "  maxO:" << maxO << ", minO:" << minO << "\n";
    cout << "  maxL:" << maxL << ", minL:" << minL << "\n";
    return result;
}

// ------------------ operator>> ------------------ 
//! read OSCAR formatted particles from a stream
istream& operator>>(istream& i, COSCARLine& line){
    i >> line.index >> line.pid 
      >> line.p[1] >> line.p[2] >> line.p[3] >> line.p[0] 
      >> line.mass 
      >> line.x[1] >> line.x[2] >> line.x[3] >> line.x[0];
    line.p[0]*=1000.0;line.p[1]*=1000.0;line.p[2]*=1000.0;line.p[3]*=1000.0;
    return i;
}
