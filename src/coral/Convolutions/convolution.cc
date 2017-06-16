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
#include <string>
#include "convolution.h"
#include "soubase.h"
#include "kernel_chooser.h"
#include "cheezyparser.h"
#include "sf.h"

#define USE_SCOTTS_KERNEL

#ifndef PI
#define PI 3.1415926535
#endif

using namespace std;
using namespace TNT;

double CConvoluterIntegrand1d::getValue(void* me, double r){
    CConvoluterIntegrand1d* metoo = (CConvoluterIntegrand1d*)me;
#ifdef USE_SCOTTS_KERNEL
    return 4.0*PI*r*r*
        (metoo->kernelPtr->GetValue(metoo->l,metoo->q,r))*
        (metoo->sourcePtr->getValue(r));
#else
    return 4.0*PI*r*r*
        Bessel::jn(metoo->l,2.0*r*(metoo->q)/197.3269602)*
        (metoo->sourcePtr->getValue(r));
#endif
}
double CConvoluterIntegrand1d::integrate(double rmin, double rmax){
    return integrator.compute(rmin,rmax);
}

CCorrFtn1dHisto convolute(const CObject1d& sourceFtn, const parameterMap& m){
    const CSourceFtnBase* sPtr = dynamic_cast<const CSourceFtnBase*>(&sourceFtn);
    string p1(""), p2("");
    if ( sPtr != NULL ) // just in case the CObject1d isn't really a source!
    { 
        p1 = sPtr->particle1;
        p2 = sPtr->particle2;
    }
    p1 = parameter::getS(m,"particle1",p1);
    p2 = parameter::getS(m,"particle2",p2);
    parameterMap convMap, corrMap;
    corrMap = parameter::getMap(m,"correlation_settings"); // manditory
    double dx = parameter::getD(corrMap,"dx");
    double xoffset = parameter::getD(corrMap,"xoffset");
    int ndata = parameter::getI(corrMap,"ndata");
    convMap = parameter::getMap(m,"convoluter_settings");// manditory
    double rmin = parameter::getD(convMap,"rmin",0.0);
    double rmax = parameter::getD(convMap,"rmax",50.0);
    cout << " _____Convolutor parameters_____"<<m<<endl<<endl;
    CKernel* kernelPtr = NULL;
    if (m.find("kernel_settings")!=m.end()) {
        kernelPtr = chooseKernel(p1,p2,parameter::getMap(m,"kernel_settings"));
    }
    else kernelPtr = chooseKernel(p1,p2,m);
    CCorrFtn1dHisto corrOut( p1, p2, false, sourceFtn.l, sourceFtn.m, true, ndata ); 
    corrOut.setFixedWidthBins(dx,xoffset);
    corrOut.covmtx_is_active=false;
    CConvoluterIntegrand1d integrand( kernelPtr, &sourceFtn );
    for (int n=0;n<corrOut.ndata;++n){
        integrand.q = corrOut.midBin(n);
        corrOut.data[n]=integrand.integrate(rmin,rmax);
        if (sourceFtn.l==0) corrOut.data[n]+=1.0;
        corrOut.uncert[n]=0.0;
    }
    return corrOut;
}
