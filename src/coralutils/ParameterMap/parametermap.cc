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
#ifndef __PARAMETERMAP__CC__
#define __PARAMETERMAP__CC__

#include "parametermap.h"
#include "cheezyparser.h"
#include "any.h"
#include "message.h"

using namespace std;

//These functions are all in the namespace parameter.
//namespace parameter {

// Getter functions with defaults
bool  parameter::getB( parameterMap m, string k, bool  d ){return m.getItem(k,d);}
int  parameter::getI( parameterMap m, string k, int  d ){return m.getItem(k,d);}
string parameter::getS( parameterMap m, string k, string d ){return m.getItem(k,d);}
double parameter::getD( parameterMap m, string k, double d ){return m.getItem(k,d);}
vector< bool >  parameter::getVB( parameterMap m, string k, const vector< bool >&  d ){return m.getItem(k,d);}
vector< int >  parameter::getVI( parameterMap m, string k, const vector< int >&  d ){return m.getItem(k,d);}
vector< double > parameter::getV( parameterMap m, string k, const vector< double >& d ){return m.getItem(k,d);}
vector< string > parameter::getVS( parameterMap m, string k, const vector< string >& d ){return m.getItem(k,d);}
vector< vector< double > > parameter::getM( parameterMap m, string k, const vector< vector< double > >& d ){return m.getItem(k,d);}

// Getter functions with no defaults
bool  parameter::getB( parameterMap m, string k ){return m.getItem<bool>(k);}
int  parameter::getI( parameterMap m, string k ){return m.getItem<int>(k);}
string parameter::getS( parameterMap m, string k ){return m.getItem<string>(k);}
double parameter::getD( parameterMap m, string k ){return m.getItem<double>(k);}
vector< bool >  parameter::getVB( parameterMap m, string k ){return m.getItem< vector< bool > >(k);}
vector< int >  parameter::getVI( parameterMap m, string k ){return m.getItem< vector< int > >(k);}
vector< double > parameter::getV( parameterMap m, string k ){return m.getItem< vector< double > >(k);}
vector< string > parameter::getVS( parameterMap m, string k ){return m.getItem< vector< string > >(k);}
vector< vector< double > > parameter::getM( parameterMap m, string k ){return m.getItem< vector< vector< double > > >(k);}
parameterMap parameter::getMap( parameterMap m, string k ){return m.getItem<parameterMap>(k);}

// Setters
void parameter::set( parameterMap& m, string k, bool  v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, int  v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, double v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, string v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, char* v ){m.setItem(k,(string)v);}

void parameter::set( parameterMap& m, string k, vector< bool >  v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, vector< int >  v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, vector< double > v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, vector< string > v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, vector< vector< double > > v ){m.setItem(k,v);}

void parameter::set( parameterMap& m, string k, parameterMap v ){m.setItem(k,v);}

//};

//-------------------------------------------------
// Reads parameterMap from a file
//-------------------------------------------------
void parameter::ReadParsFromFile(parameterMap& apars, string inFileName){
	ifstream inMapFile(inFileName.c_str());
	if (!inMapFile.good())
		throw MESSAGE<< "Cannot open file "<<inFileName<<ENDM_SEVERE;
	inMapFile>>apars;
}

void parameter::ReadParsFromFile(parameterMap& apars, char *inFileName){
	cout << "Reading parameters from file " << inFileName << endl;
	ifstream inMapFile(inFileName);
	if (!inMapFile.good())
		throw MESSAGE<< "Cannot open file "<<inFileName<<ENDM_SEVERE;
	inMapFile>>apars;
}
//-------------------------------------------------
// Writes parameterMap to a file
//-------------------------------------------------
void parameter::WriteParsToFile(parameterMap& apars, string outFileName){
	ofstream outMapFile(outFileName.c_str());
	if (!outMapFile.good())
		throw MESSAGE<< "Cannot open file "<<outFileName<<ENDM_SEVERE;
	outMapFile<<apars;
}

void parameter::PrintPars(parameterMap& pars){
	cout << pars << endl;
	/*
	map<string,boost::any>::iterator itr;
	for(itr=pars.begin(); itr !=pars.end(); ++itr){
		cout << itr->first << " = " << itr->second << endl;
	}*/
}

#endif
