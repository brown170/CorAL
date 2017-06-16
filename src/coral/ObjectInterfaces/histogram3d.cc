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
#include "histogram3d.h"
#include "message.h"

// Constructors & Destructors
CHistogram3d::CHistogram3d(int Nx,int Ny,int Nz,double Dx,double Dy,double Dz,double xO, double yO, double zO):
    CDataSet(Nx*Ny*Nz), CObject3d(), nx(Nx), ny(Ny), nz(Nz), dx(Dx), dy(Dy), dz(Dz),
    xoffset(xO), yoffset(yO), zoffset(zO) {
        ixzero=findBin(dx/2.,dx,xoffset);
        iyzero=findBin(dy/2.,dy,yoffset);
        izzero=findBin(dz/2.,dz,zoffset);
    }

CHistogram3d::CHistogram3d(const CHistogram3d& a): 
    CDataSet(a), CObject3d(a),
    nx(a.nx),ny(a.ny),nz(a.nz),
    dx(a.dx),dy(a.dy),dz(a.dz),
    xoffset(a.xoffset),yoffset(a.yoffset),zoffset(a.zoffset),
    ixzero(a.ixzero),iyzero(a.iyzero),izzero(a.izzero){}

// Copy Functions
void CHistogram3d::CopyState(const CHistogram3d& a){
    CDataSet::CopyState(a);
    CObject3d::CopyState(a);
    nx=a.nx;ny=a.ny;nz=a.nz;
    dx=a.dx;dy=a.dy;dz=a.dz;
    xoffset=a.xoffset;yoffset=a.yoffset;zoffset=a.zoffset;
    ixzero=a.ixzero;iyzero=a.iyzero;izzero=a.izzero;
}
    
// Bin index finding
bool CHistogram3d::inThisBin(int ix, int iy, int iz, double x, double y, double z) const{
    return whatBin(x,y,z)==whatIndex(ix,iy,iz);
}

int CHistogram3d::whatBin(double x, double y, double z) const {
    int i=findBin(x,dx,xoffset);
    int j=findBin(y,dy,yoffset);
    int k=findBin(z,dz,zoffset);
    if (i>=0 && i<nx && j>=0 && j<ny && k>=0 && k<nz) return whatIndex(i,j,k);
    else return -1;
}

double CHistogram3d::getValueCart(double x, double y, double z) const {
    int iBin = whatBin(x,y,z);
    if ( iBin < 0 ) return 0.0;
    return data[iBin];
}

double CHistogram3d::getErrorCart(double x, double y, double z) const {
    int iBin = whatBin(x,y,z);
    if ( iBin < 0 ) return 0.0;
    return uncert[iBin];
}

// i/o
bool CHistogram3d::Read(const parameterMap& m){
    CObject3d::Read(m);
    if (m.find("crab_style_datablock")!=m.end()){
        throw MESSAGE<<"write me! crab_style_datablock"<<ENDM_FATAL;
        return false;
    } else { // do it easy way
        nx=parameter::getI(m,"nx",nx);
        ny=parameter::getI(m,"ny",ny);
        nz=parameter::getI(m,"nz",nz);
        redim( nx*ny*nz );
        CDataSet::Read(m);
        if (ndata != nx*ny*nz) {
            MESSAGE << "ndata = "<<ndata<<" while nx*ny*nz = "<<nx*ny*nz<<", so we are trashing the data in the parameterMap"<<ENDM_WARN;
            redim( nx*ny*nz );
        }
        dx=parameter::getD(m,"dx",dx);
        dy=parameter::getD(m,"dy",dy);
        dz=parameter::getD(m,"dz",dz);
        xoffset=parameter::getD(m,"xoffset",xoffset);
        yoffset=parameter::getD(m,"yoffset",yoffset);
        zoffset=parameter::getD(m,"zoffset",zoffset);
        ixzero=findBin(dx/2.,dx,xoffset);
        iyzero=findBin(dy/2.,dy,yoffset);
        izzero=findBin(dz/2.,dz,zoffset);
    }
    return true;
}

bool CHistogram3d::Write(parameterMap& m){
    CObject3d::Write(m);
    parameter::set(m,"nx",nx);
    parameter::set(m,"ny",ny);
    parameter::set(m,"nz",nz);
    parameter::set(m,"dx",dx);
    parameter::set(m,"dy",dy);
    parameter::set(m,"dz",dz);
    parameter::set(m,"xoffset",xoffset);
    parameter::set(m,"yoffset",yoffset);
    parameter::set(m,"zoffset",zoffset);
    return CDataSet::Write(m);
}
    
