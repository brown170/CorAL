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
#ifndef XGRAPH_INCLUDE_H
#define XGRAPH_INCLUDE_H
#include <unistd.h>
#include <iostream>
#include <X11/Xlib.h>
#include <cmath>
#include <cstring>
#include <cstdlib>
//#include <dlfcn.h>
#include <unistd.h>
//#include <ltdl.h>

using namespace std;

class CAxesInfo{
 public:
  double xmin,ymin,xmax,ymax;
  double xscale,yscale;
  double xintercept; //where y-axis intercepts x axis
  double yintercept; //where x-axis intercepts y axis
  int nxtics;
  int nytics;
  CAxesInfo();
};

class CXGraph{
 public:
  Display *theDisplay;
  Window theWindow;
  Screen *theScreen;
  Font theFont;
  GC theGC;
	int window_horizsize,window_vertsize,ixposition,iyposition;
  Colormap screen_colormap;
  XColor red, brown, blue, yellow, green, orange, black, cyan, violet, lightblue, navy, pink,white;
  Status rc;
  CAxesInfo axesinfo;
  void setcolor(string color);
  void setaxes(double xmin,double ymin,double xmax,double ymax);
  void drawtext(char *string,double x,double y);
  void drawline(double x1,double y1,double x2,double y2);
  void drawcircle(double x,double y,double size);
  void drawsquare(double x,double y,double size);
  void drawdiamond(double x,double y,double size);
  void drawuptriangle(double x,double y,double size);
  void drawdowntriangle(double x,double y,double size);
  void drawpoint(double x,double y);
  void drawarrow(double x1,double y1,double x2,double y2,double headsize);
  void drawaxes();
  void closedisplay();
  void getij(double x,double y,int *i,int *j);
  void plotline(double *x,double *y,int npts);
  void plotpoints(double *x,double *y,int npts);
	void clear();
	void flush();
  CXGraph(int horizsize,int vertsize,int ixposition,int iyposition);
	~CXGraph();
};
#endif
