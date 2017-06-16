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
//#include <stream>
#include <cstdlib>
#include <cmath>
//#include <dlfcn.h>
#include <X11/Xlib.h>
#include "xgraph.h"
#include <cstring>

const double PI=3.14159265358979323844;

int main(){
  double xmin,ymin,xmax,ymax;
  double x,y;
  char dummy[100];
  int window_xoff=20,window_yoff=100,window_width=500,window_height=500;
  double x1,y1,x2,y2;
  double symbolsize=0.02;
  char title[20];
  int i,npts;

  xmin=0.0;
  ymin=-1.1;
  xmax=12.0;
  ymax=1.1;
  //printf("Enter width, height, xoff, yoff :");
  //scanf("%d %d %d %d",&window_width,&window_height,&window_xoff,&window_yoff);
  CXGraph xgraph(window_width,window_height,window_xoff,window_yoff);
  window_xoff+=window_width+20;
  CXGraph ygraph(window_width,window_height,window_xoff,window_yoff);

  xgraph.setaxes(xmin,ymin,xmax,ymax);
  xgraph.drawaxes();

  x=0.5*(xmax+xmin);
  y=0.5*(ymin+ymax);
  xgraph.drawtext(title,x,y);
  xgraph.setcolor("cyan");

  
  ygraph.setaxes(xmin,ymin,xmax,ymax);
  ygraph.drawaxes();
  x1=xmin+0.2*(xmax-xmin);
  y1=ymin+0.2*(ymax-ymin);
  x2=xmin+0.4*(xmax-xmin);
  y2=ymin+0.6*(ymax-ymin);
  ygraph.drawarrow(x1,y1,x2,y2,.03);

  npts=50;
  for(i=0;i<=npts;i++){
    x=double(i)*xmax/double(npts);
    y=sin(x);
    xgraph.drawpoint(x,y);
    xgraph.drawcircle(x,y,symbolsize);
    xgraph.drawsquare(x,y,symbolsize);
    xgraph.drawdiamond(x,y,symbolsize);
    xgraph.drawuptriangle(x,y,symbolsize);
    xgraph.drawdowntriangle(x,y,symbolsize);
  }


  //pause();
  printf("Enter anything :");
  scanf("%s",&dummy);
  xgraph.closedisplay();
  ygraph.closedisplay();
  return 0;
}

#include "xgraph.cc"
