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
#ifndef XGRAPH_INCLUDE_CC
#define XGRAPH_INCLUDE_CC
#include "xgraph.h"

CXGraph::CXGraph(int horizsize,int vertsize,int ixposition_in,int iyposition_in){
	window_horizsize=horizsize;
	window_vertsize=vertsize;
	ixposition=ixposition_in;
	iyposition=iyposition_in;
	theDisplay = XOpenDisplay(NULL);
	XSynchronize(theDisplay, True);
	theScreen = DefaultScreenOfDisplay(theDisplay);
	theWindow = XCreateSimpleWindow(theDisplay, RootWindowOfScreen(theScreen), 
		ixposition,iyposition, 
		window_horizsize, window_vertsize, 0, 
		BlackPixelOfScreen(theScreen), 
		WhitePixelOfScreen(theScreen));
	theGC = XCreateGC(theDisplay, theWindow, 0L, NULL);
	XSetForeground(theDisplay, theGC, BlackPixelOfScreen(theScreen));
	XMapWindow(theDisplay,theWindow);
	XMoveWindow(theDisplay,theWindow,ixposition+50,iyposition+50);
	XSetLineAttributes(theDisplay,theGC,2,
		LineSolid,CapRound,JoinRound);
	screen_colormap = DefaultColormap(theDisplay, DefaultScreen(theDisplay));

	rc = XAllocNamedColor(theDisplay, screen_colormap, "red", &red, &red);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'red' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "brown", &brown, &brown);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'brown' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "blue", &blue, &blue);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'blue' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "yellow", &yellow, &yellow);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'yellow' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "green", &green, &green);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'green' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "orange", &orange, &orange);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'orange' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "black", &black, &black);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'black' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "cyan", &cyan, &cyan);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'cyan' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "violet", &violet, &violet);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'violet' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "lightblue", &lightblue, &lightblue);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'lightblue' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "navy", &navy, &navy);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'navy' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "pink", &pink, &pink);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'pink' color.\n");
		exit(1);
	}
	rc = XAllocNamedColor(theDisplay, screen_colormap, "white", &white, &white);
	if (rc == 0) {
		fprintf(stderr, "XAllocNamedColor - failed to allocated 'white' color.\n");
		exit(1);
	}

	//theFont=XLoadFont(theDisplay,"9x15");
}

CXGraph::~CXGraph(){
	XFlush(theDisplay);
}
void CXGraph::flush(){
	XFlush(theDisplay);
	XCloseDisplay(theDisplay);
}

void CXGraph::setcolor(string color){
	if(color=="red"){
		XSetForeground(theDisplay, theGC, red.pixel);
	}
	if(color=="brown"){
		XSetForeground(theDisplay, theGC, brown.pixel);
	}
	if(color=="yellow"){
		XSetForeground(theDisplay, theGC, yellow.pixel);
	}
	if(color=="blue"){
		XSetForeground(theDisplay, theGC, blue.pixel);
	}
	if(color=="green"){
		XSetForeground(theDisplay, theGC, green.pixel);
	}
	if(color=="orange"){
		XSetForeground(theDisplay, theGC, orange.pixel);
	}
	if(color=="black"){
		XSetForeground(theDisplay, theGC, black.pixel);
	}
	if(color=="cyan"){
		XSetForeground(theDisplay, theGC, cyan.pixel);
	}
	if(color=="violet"){
		XSetForeground(theDisplay, theGC, violet.pixel);
	}
	if(color=="lightblue"){
		XSetForeground(theDisplay, theGC, lightblue.pixel);
	}
	if(color=="navy"){
		XSetForeground(theDisplay, theGC, navy.pixel);
	}
	if(color=="pink"){
		XSetForeground(theDisplay, theGC, pink.pixel);
	}
}

void CXGraph::drawtext(char *charstring,double x,double y){
	int i,j;
	getij(x,y,&i,&j);
	i=i-strlen(charstring)*3;
	XDrawString(theDisplay,theWindow,theGC,i,j,charstring,strlen(charstring));
}

CAxesInfo::CAxesInfo(){
	xintercept=yintercept=0.0;
	nxtics=nytics=4;
}

void CXGraph::drawline(double x1,double y1,double x2,double y2){
	int i1,j1,i2,j2;
	getij(x1,y1,&i1,&j1);
	getij(x2,y2,&i2,&j2);
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}

void CXGraph::drawarrow(double x1,double y1,double x2,double y2,
double headsize){
	int i1,j1,i2,j2,ihead,iperp,jperp,ipar,jpar;
	double norm;
	getij(x1,y1,&i1,&j1);
	getij(x2,y2,&i2,&j2);
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	ihead=int(headsize*0.8*window_horizsize);
	norm=sqrt(double((i1-i2)*(i1-i2)+(j1-j2)*(j1-j2)));
	ipar=int(ihead*(i1-i2)/norm);
	jpar=int(ihead*(j1-j2)/norm);
	iperp=jpar;
	jperp=-ipar;
	i1=i2+ipar+iperp;
	j1=j2+jpar+jperp;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i1=i1-2*iperp;
	j1=j1-2*jperp;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}

void CXGraph::drawpoint(double x,double y){
	int i,j;
	getij(x,y,&i,&j);
	XDrawPoint(theDisplay,theWindow,theGC,i,j);
}

void CXGraph::drawcircle(double x,double y,double size){
	int i,j,ir;
	const double PI=4.0*atan(1.0);
	getij(x,y,&i,&j);
	ir=1+int(2.0*size*0.8*window_horizsize/sqrt(PI));
	XDrawArc(theDisplay,theWindow,theGC,i-ir,j-ir,2*ir,2*ir,0,23040);
}

void CXGraph::drawsquare(double x,double y,double size){
	int i1,j1,i2,j2,i,j,h;
	getij(x,y,&i,&j);
	h=int(size*0.8*window_horizsize);
	i1=i-h;
	j1=j-h;
	i2=i1+2*h;
	j2=j1;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i1;
	j2=j1+2*h;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i1=i2+2*h;
	j1=j2;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i1;
	j2=j1-2*h;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}

void CXGraph::drawdiamond(double x,double y,double size){
	int i1,j1,i2,j2,i,j,h;
	h=int(size*0.8*window_horizsize*sqrt(2.0));
	getij(x,y,&i,&j);
	i1=i;
	j1=j+h;
	i2=i-h;
	j2=j;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i+h;
	j2=j;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i1=i;
	j1=j-h;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i-h;
	j2=j;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}

void CXGraph::drawuptriangle(double x,double y,double size){
	int i1,j1,i2,j2,i,j,h;
	h=int(size*0.8*window_horizsize*sqrt(2.0));
	getij(x,y,&i,&j);
	i1=i;
	j1=j+int(h*(sqrt(3.0)-1.0/sqrt(3.0)));
	i2=i-h;
	j2=j-int(h/sqrt(3.0));
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i+h;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i1=i-h;
	j1=j-int(h/sqrt(3.0));
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}

void CXGraph::drawdowntriangle(double x,double y,double size){
	int i1,j1,i2,j2,i,j,h;
	getij(x,y,&i,&j);
	h=int(size*0.8*window_horizsize*sqrt(2.0));
	i1=i;
	j1=j-int(h*(sqrt(3.0)-1.0/sqrt(3.0)));
	i2=i-h;
	j2=j+int(h/sqrt(3.0));
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i2=i+h;
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
	i1=i-h;
	j1=j+int(h/sqrt(3.0));
	XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
}


void CXGraph::setaxes(double xmin,double ymin,double xmax,double ymax){
	axesinfo.xmin=xmin;
	axesinfo.ymin=ymin;
	axesinfo.xmax=xmax;
	axesinfo.ymax=ymax;
	axesinfo.xscale =xmax-xmin;
	axesinfo.yscale=ymax-ymin;
}

void CXGraph::drawaxes(){
	int itic;
	double x1,x2,y1,y2;
	int i1,j1,i2,j2;
	int arrowsize,ticsize;
	double xtext,ytext;
	char charstring[40];
	ticsize=int(0.02*0.8*window_horizsize);
	arrowsize=ticsize;

	x1=axesinfo.xintercept;
	y1=axesinfo.ymin;
	x2=x1;
	y2=axesinfo.ymax+0.06*axesinfo.yscale;
	drawarrow(x1,y1,x2,y2,0.02);
	for(itic=0;itic<=axesinfo.nytics;itic++){
		y1=axesinfo.ymin+itic*axesinfo.yscale/double(axesinfo.nytics);
		x1=axesinfo.xintercept;
		getij(x1,y1,&i1,&j1);
		i1=i1-ticsize;
		i2=i1+2*ticsize;
		j2=j1;
		XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
		xtext=x1-30*(axesinfo.xscale)/(0.8*window_vertsize);
		ytext=y1-4*(axesinfo.yscale)/(0.8*window_vertsize);
		sprintf(charstring,"%g",y1);
		drawtext(charstring,xtext,ytext);
	}  

	y1=axesinfo.yintercept;
	x1=axesinfo.xmin;
	y2=y1;
	x2=axesinfo.xmax+0.06*axesinfo.xscale;
	drawarrow(x1,y1,x2,y2,.02);
	for(itic=0;itic<=axesinfo.nxtics;itic++){
		x1=axesinfo.xmin+itic*axesinfo.xscale/double(axesinfo.nxtics);
		y1=axesinfo.yintercept;
		getij(x1,y1,&i1,&j1);
		j1=j1-ticsize;
		j2=j1+2*ticsize;
		i2=i1;
		XDrawLine(theDisplay,theWindow,theGC,i1,j1,i2,j2);
		xtext=x1;
		ytext=y1-20*(axesinfo.yscale)/(0.8*window_vertsize);
		sprintf(charstring,"%g",xtext);
		drawtext(charstring,xtext,ytext);
	}
}

void CXGraph::getij(double x,double y,int *i,int *j){
	*i=int(0.15*window_horizsize
		+0.8*window_horizsize*(x-axesinfo.xmin)
		/(axesinfo.xmax-axesinfo.xmin));
	*j=int(0.9*window_vertsize
		-0.8*window_vertsize*(y-axesinfo.ymin)
		/(axesinfo.ymax-axesinfo.ymin));
}

void CXGraph::closedisplay(){
	XCloseDisplay(theDisplay);
}

void CXGraph::plotline(double *x,double *y,int npts){
	double xx,yy,oldx,oldy;
	double xmin,ymin,xmax,ymax;
	int i;
	xmin=axesinfo.xmin;
	ymin=axesinfo.ymin;
	xmax=axesinfo.xmax;
	ymax=axesinfo.ymax;

	for(i=1;i<npts;i++){
		xx=x[i]; yy=y[i];
		oldx=x[i-1]; oldy=y[i-1];
		if(xx<=xmax && xx>=xmin && yy<=ymax && yy>=ymin
		&& oldx<=xmax && oldx>=xmin && oldy <=ymax && oldy>=ymin){
			drawline(oldx,oldy,xx,yy);
		}
	}
}

void CXGraph::plotpoints(double *x,double *y,int npts){
	double xx,yy,size;
	double xmin,ymin,xmax,ymax;
	int i;
	xmin=axesinfo.xmin;
	ymin=axesinfo.ymin;
	xmax=axesinfo.xmax;
	ymax=axesinfo.ymax;
	size=0.005;

	for(i=0;i<npts;i++){
		xx=x[i]; yy=y[i];
		if(xx<=xmax && xx>=xmin && yy<=ymax && yy>=ymin){
			drawsquare(xx,yy,size);
		}
	}
}

void CXGraph::clear(){
  XSetForeground(theDisplay,theGC,white.pixel);
	XFillRectangle(theDisplay,theWindow,theGC,0,0,window_horizsize,window_vertsize);
  XSetForeground(theDisplay,theGC,black.pixel);
}
#endif
