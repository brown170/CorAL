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
#ifndef HISTOGRAM3D_H
#define HISTOGRAM3D_H

#include <cmath>
#include "parametermap.h"
#include "dataset.h"
#include "objects3d.h"
#include "tnt_array1d.h"
#include "utils.h"

//----------------------------------------------
// Simple 3d histogram class, only fixed width bins supported
//----------------------------------------------
class CHistogram3d: public CDataSet, public CObject3d {

public:

    int nx,ny,nz;
    double dx,dy,dz;
    double xoffset,yoffset,zoffset;
    int ixzero,iyzero,izzero;
  
    // Constructors & Destructors
    CHistogram3d(int Nx=0,int Ny=0,int Nz=0,
        double Dx=1.,double Dy=1.,double Dz=1.,
        double xO=0.5, double yO=0.5, double zO=0.5);
    CHistogram3d(const CHistogram3d& a);

    // i/o
    bool Read(const parameterMap& m);
    bool Write(parameterMap& m);
    
    // Copy Functions
    void CopyState(const CHistogram3d& a);
    
    // Main interface to source value 
    virtual double getValueCart(double x, double y, double z) const;
    virtual double getValueSphr(double r, double theta, double phi) const 
        {double rho=r*sin(theta);return getValueCart(rho*cos(phi),rho*sin(phi),r*cos(theta));}
    virtual double getErrorCart(double x, double y, double z) const;
    virtual double getErrorSphr(double r, double theta, double phi) const 
        {double rho=r*sin(theta);return getErrorCart(rho*cos(phi),rho*sin(phi),r*cos(theta));}

    // Bin index finding
    int findBin(double x, double _dx, double _xoffset) const{return iround<double>((x-_xoffset)/_dx);}
    int whatIndex(int ix, int iy, int iz) const{return iz*ny*nx+iy*nx+ix;}
    bool inThisBin(int ix, int iy, int iz, double x, double y, double z) const;
    int whatBin(double x, double y, double z) const;

    // Centers of bins
    double midBinX(int ix) const{return ix*dx+xoffset;}
    double midBinY(int iy) const{return iy*dy+yoffset;}
    double midBinZ(int iz) const{return iz*dz+zoffset;}
    Array1D<double> binCenter(int ix,int iy,int iz) const
        {Array1D<double> r(3); r[0]=midBinX(ix); r[1]=midBinY(iy); r[2]=midBinZ(iz); return r;}
    Array1D<double> binCenter(int i) const
        {Array1D<double> r(3); r[0]=midBinX((i%(ny*nx))%nx); r[1]=midBinY((i%(ny*nx))/nx); r[2]=midBinZ(i/(ny*nx)); return r;}

    // Left edges of bins (to get right edge, put in i+1)
    double leftBinEdgeX(int ix) const{return midBinX(ix)-dx/2.;}
    double leftBinEdgeY(int iy) const{return midBinY(iy)-dy/2.;}
    double leftBinEdgeZ(int iz) const{return midBinZ(iz)-dz/2.;}
    // Right edges of bins (to get left edge, put in i-1)
    double rightBinEdgeX(int ix) const{return midBinX(ix)+dx/2.;}
    double rightBinEdgeY(int iy) const{return midBinY(iy)+dy/2.;}
    double rightBinEdgeZ(int iz) const{return midBinZ(iz)+dz/2.;}

    // Find index of reflected bins
    int flipBin(int i, int _izero) const{int ii=2*_izero-1-i;if(ii<0){return i;}else{return ii;}} 
    int flipBinX(int ix) const{return flipBin(ix,ixzero);}
    int flipBinY(int iy) const{return flipBin(iy,iyzero);}
    int flipBinZ(int iz) const{return flipBin(iz,izzero);}

    // Misc. interfaces
    double binVolume(void) const{return dx*dy*dz;}
};

#endif
