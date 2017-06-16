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
#include "histogram1d.h"
#include "message.h"
#include "linalg.h"

void CHistogram1d::setFixedWidthBins(double dx, double xoffset){
    Array1D<double> tmp(ndata+1,0.);
    for (int i=0;i<ndata+1;++i) {tmp[i]=dx*i+xoffset-dx/2.;}
    knots=tmp;
    fixed_width_bins=true;
    xmin = getLeftSupport( 0 );
    xmax = getRightSupport( ndata - 1 );
}
    
bool CHistogram1d::inThisBin(int i, double xx) const
    {return ((xx<knots[i+1])&&(xx>=knots[i]));}
    
//! Computes which bin x lies in
//! The graceful flag determines what happens if x is off grid:
//!   If graceful==true, give the nearest on-grid bin
//!   If graceful==false, give the nearest off-grid "bin".  You must check for this as
//!     it mostlikely will cause an assertion error if bounds checking is on or a segfault
//!     if bounds checking is off. 
int CHistogram1d::whatBin(double x, bool graceful) const{
    for (int i=0;i<ndata;++i){if (inThisBin(i,x)) return i;} 
    if (x<knots[0]) {
        MESSAGE << "x = "<<x<<" is off scale, "<< "min x="<<knots[0];
        if (graceful) MESSAGE<<ENDM_WARN; else 
            MESSAGE<<ENDM_INFO;
        if (graceful) return 0; else return -1;
    } else {
        MESSAGE << "x = "<<x<<" is off scale, "<< "max x="<<knots[ndata];
        //if (graceful) MESSAGE<<ENDM_WARN; else 
            MESSAGE<<ENDM_INFO;
        if (graceful) return ndata-1; else return ndata;
    }
}
    
double CHistogram1d::basisFunction(double x, int i, int jderiv) const{
    switch (jderiv){
        case 0:  return double(inThisBin(i,x));
        case -1: return x*double(inThisBin(i,x));
        default: return 0.0;
    }
    return 0.0; // should never get here!
}
    
// i/o
bool CHistogram1d::Read(const parameterMap& m){
    bool found_bins=false, found_data=false;
    // some defaults
    vector< double > empty_vec(0);
    vector< vector< double > > empty_mat(0);
    // read base-class stuff
    CObject1d::Read(m);
    // get the binning data if not in datablock
    fixed_width_bins=parameter::getB(m,"fixed_width_bins",fixed_width_bins);
    if (fixed_width_bins){
        found_bins=true;
        CDataSet::Read(m);
        if ( m.hasKey("dx") && m.hasKey("xoffset") ) setFixedWidthBins(parameter::getD(m,"dx"),parameter::getD(m,"xoffset"));
        else if ( m.hasKey("xmin") && m.hasKey("xmax") ) {
            double binWidth = ( parameter::getD(m,"xmax") - parameter::getD(m,"xmin") )/ndata;
            setFixedWidthBins(binWidth,binWidth/2.0+m.hasKey("xmin"));
        }
        else throw MESSAGE << "Cannot find dx or xoffset keys or xmax and xmin keys"<<ENDM_FATAL;
    } else if(m.find("left_bin_edges")!=m.end()){
        found_bins=true;
        CDataSet::Read(m);
        knots = stl2tntVec(parameter::getV(m,"knots",empty_vec));
        xmin = getLeftSupport( 0 );
        xmax = getRightSupport( ndata - 1 );
    } 
    // read the data itself using the CDataSet base class if possible
    found_data=CDataSet::Read(m);
    // for datablock types not supported by base class, we do it ourself
    if(m.find("xy_datablock")!=m.end()) {
        found_bins=true;
        found_data=true;
        vector< vector<double> > tmp;
        tmp=parameter::getM(m,"xy_datablock",empty_mat);
        ndata = tmp.size();
        Array1D<double> xtmp(ndata+1,0.);
        Array1D<double> ytmp(ndata,0.);
        for (int i=0;i<ndata;++i){
            xtmp[i+1]=2.*tmp[0][i]-xtmp[i];
            ytmp[i]=tmp[1][i];
        }
        data = ytmp;
        if (!fixed_width_bins) knots = xtmp;
        xmin = getLeftSupport( 0 );
        xmax = getRightSupport( ndata - 1 );
    } else if(m.find("xydy_datablock")!=m.end()) {
        found_bins=true;
        found_data=true;
        vector< vector<double> > tmp;
        tmp=parameter::getM(m,"xydy_datablock",empty_mat);
        ndata = tmp.size();
        Array1D<double> xtmp(ndata+1,0.);
        Array1D<double> ytmp(ndata,0.);
        Array1D<double> dytmp(ndata,0.);
        for (int i=0;i<ndata;++i){
            xtmp[i+1]=2.*tmp[0][i]-xtmp[i];                      
            ytmp[i]=tmp[1][i];
            dytmp[i]=tmp[2][i];
        }
        data = ytmp;
        uncert=dytmp;
        if (!fixed_width_bins) knots = xtmp;
        xmin = getLeftSupport( 0 );
        xmax = getRightSupport( ndata - 1 );
    }
    if (!found_bins) {
        MESSAGE << "Bin edges not found, setting to zero!"<<ENDM_INFO;
        Array1D<double> xtmp(ndata+1,0.);
        knots=xtmp;
    }
    if (!found_data) {
        MESSAGE << "Data set not found, setting to zero!" <<ENDM_INFO;
        Array1D<double> ytmp(ndata,0.), dytmp(ndata,0.);
        data=ytmp;
        uncert=dytmp;
    }
    return true;
}

bool CHistogram1d::Write(parameterMap& m){
    CObject1d::Write(m);
    if (fixed_width_bins) {
        parameter::set(m,"fixed_width_bins",true);
        parameter::set(m,"dx",(knots[1]-knots[0]));
        parameter::set(m,"xoffset",(knots[1]-knots[0])/2.);
    } else {
        parameter::set(m,"left_bin_edges",tnt2stlVec(knots));
    }
    return CDataSet::Write(m);
}

bool CHistogram1d::setDefaultKnots( void ){
    double dx=(xmax-xmin)/(double)ndata; 
    setFixedWidthBins(dx,dx/2.0); 
    return true;
}
