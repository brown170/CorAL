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
#include <cmath>
#include "dataset.h"
#include "linalg.h"
#include "message.h"

/// read to parameter map
bool CDataSet::Read(const parameterMap& com){

    bool found_data=false;

    // figure out the number of data points
    ndata=parameter::getI(com,"ndata",data.dim());

    vector< vector < double > > empty_mat(0);
    vector< double > empty_vec(0);

    // load any covariance matrix
    if (com.find("covmtx_datablock")!=com.end()) {
        MESSAGE << "loading covmtx_datablock"<<ENDM_INFO;
        covmtx=stl2tntMat(parameter::getM(com,"covmtx_datablock",empty_mat));
        covmtx_is_active=true;
        syncUncert();
    } else {
        covmtx_is_active=false;
    }


    // check for and load any data
    if ((com.find("ydy_datablock")==com.end())&&(com.find("y_datablock")==com.end())){
        MESSAGE << "Data set not found, setting to zero!" <<ENDM_INFO;
        Array1D<double> ytmp(ndata,0.), dytmp(ndata,0.);
        data=ytmp;
        uncert=dytmp;
        return false;
    }

    //OK two possiblities from here...
    
    // first, just data in map
    if (com.find("y_datablock")!=com.end()) {
        found_data=true;
        MESSAGE << "Loading y_datablock"<<ENDM_INFO;
        data = stl2tntVec(parameter::getV(com,"y_datablock",empty_vec));
        if (!covmtx_is_active) {Array1D<double> tmp(data.dim(),0.); uncert=tmp;} //otw. set already
    }
    
    // second, data and uncertanties in map
    if (com.find("ydy_datablock")!=com.end()) {
        found_data=true;
        MESSAGE << "Loading ydy_datablock"<<ENDM_INFO;
        vector< vector<double> > tmp;
        tmp = parameter::getM(com,"ydy_datablock",empty_mat);
        Array1D<double> y(tmp.size(),0.), dy(tmp.size(),0.);
        if (covmtx_is_active) {
            for (int i=0;i<y.dim();++i) {y[i]=tmp[i][0];}
            syncUncert();
            data = y;
        } else {
            for (int i=0;i<y.dim();++i) {y[i]=tmp[i][0];dy[i]=abs(tmp[i][1]);}
            data = y;
            uncert = dy;
        }
    }

    // finally make sure all dimensions match, otw. we're hosed
    if (!checkArrayDimensions()) exit(-1);
    return found_data;
}

/// write to parameter map
bool CDataSet::Write(parameterMap& com){
    parameter::set(com,"ndata",ndata);
    if (covmtx_is_active){
        parameter::set(com,"covmtx_datablock",tnt2stlMat(covmtx));
        parameter::set(com,"y_datablock",tnt2stlVec(data));
    } else {
        Array2D<double> tmp(ndata,2,0.);
        for (int i=0;i<ndata;++i){tmp[i][0]=data[i];tmp[i][1]=uncert[i];}
        parameter::set(com,"ydy_datablock",tnt2stlMat(tmp));
    }
    return true;
}


bool CDataSet::checkArrayDimensions(void){
    bool result=true;
    if (data.dim() != ndata) {
        MESSAGE << "data array dim ("<<data.dim1()<<") != ndata ("<<ndata<<")!"<<ENDM_FATAL;
        result=false;
    }
    if (uncert.dim() != ndata) {
        MESSAGE << "uncert array dim ("<<uncert.dim1()<<") != ndata("<<ndata<<")!"<<ENDM_FATAL;
        result=false;
    }
    if (covmtx_is_active && (covmtx.dim1() != ndata)) {
        MESSAGE << "covmtx array dim1 ("<<covmtx.dim1()<<") != ndata("<<ndata<<")!"<<ENDM_FATAL;
        result=false;
    }
    if (covmtx_is_active && (covmtx.dim2() != ndata)) {
        MESSAGE << "covmtx array dim2 ("<<covmtx.dim2()<<") != ndata("<<ndata<<")!"<<ENDM_FATAL;
        result=false;
    }
    return result;
}

/// Redimensions all internal arrays and zeros them (use with caution!)
bool CDataSet::redim(int N){
    ndata=N;
    Array1D<double> y(N,0.),dy(N,0.);
    data=y;uncert=dy;
    if (covmtx_is_active) {
        Array2D<double> cm(N,N,0.);
        covmtx=cm;
    } else {
        Array2D<double> cm(1,1,0.);
        covmtx=cm;
    }
    return true;
}

/// Builds the error vector from the covariance matrix
void CDataSet::syncUncert(void){
    Array1D<double> u(ndata);
    for (int i=0;i<ndata;i++) u[i]=sqrt(abs(covmtx[i][i]));
    uncert=u;
}

/// Builds the covariance matrix from the error vector
void CDataSet::syncCovMtx(void){
    Array2D<double> cov(ndata,ndata,0.);
    for (int i=0;i<ndata;++i){
        cov[i][i]=uncert[i]*uncert[i];
    }
    *const_cast<bool*>(&covmtx_is_active)=true;
    covmtx=cov;
}


