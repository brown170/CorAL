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
#ifndef __NEW_PARAMETERMAP_H
#define __NEW_PARAMETERMAP_H

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
#include <iostream>
#include "message.h"
#include "convert.h"
#include "any.h"

using namespace std;

//---------------------------------------------------------------
//This code helps one parse the map that contains configuration 
// parameters.
//
//Example:
// 
// MyCoolClass(parameterMap m){  //m has the parameters and is passed in
// int importantParameter = parameter::getI(m,"nameOfParameter",-1);
//  ...
//
//  If "nameOfParameter is not a key in the map, -1 is returned since that
// is the 3rd argment of the getI function.
//
//MH 22 jun04
//---------------------------------------------------------------

//This code only works with a map of the type below.  The type def is
//to make it easy to remember.
class  parameterMap: public map< string, boost::any > {
	
	public: 
	
	template< typename ValueType > 
	ValueType getItem( const string key ){
	    if ( !hasKey( key ) ) throw MESSAGE << "Key " << key << " not found in parameterMap" << ENDM_SEVERE;
		ValueType answer;
		try {
                answer = boost::any_cast< ValueType >( find( key )->second );
		}
		catch ( boost::bad_any_cast ex ){ 
		        cout << endl;
                throw MESSAGE << ex.what() << " when converting value for key = " << key << ENDM_SEVERE;
		}
		return answer; 
	}
	
	template< typename ValueType >
	ValueType getItem( const char* key ){ return getItem< ValueType >( string( key ) ); }
	
	template< typename ValueType > 
	ValueType getItem( const string key, const ValueType& def ){
	  if ( find( key ) != end() ){ 
	    return getItem< ValueType >( key ); 
	  } else { 
	    #ifdef WARN_ON_DEFAULT
	    stringstream mymessage;
	    mymessage << "Default value of " << def << " will be used for key " << key << " for ALL reads from this map!";
	    MESSAGE << mymessage.str() << ENDM_WARN;
	    setItem( key, def );
	    #endif
	    return def;
	  }
	}
	
	template< typename ValueType > 
        ValueType getItem( const char* key, const ValueType& def ){ return getItem( string( key ), def ); }
	
	template< typename ValueType > 
	void setItem( const string key, const ValueType& val ){
		if (find(key)!=end()) erase(key); 
		insert(make_pair(key,val));
	}
	 
	template< typename ValueType > 
	void setItem( const char* key, const ValueType& val ){ return setItem( string( key ), val ); }
	
	bool hasKey( const string key ) const{ return find(key)!=end(); }

	bool hasKey( const char* key ) const{ return hasKey( string( key ) ); }
	
};

//These functions are all in the namespace parameter.
namespace parameter{
	
	// Getter functions with defaults
	bool   getB( parameterMap m, string k, bool   d );
	int    getI( parameterMap m, string k, int    d );
	string getS( parameterMap m, string k, string d );
	double getD( parameterMap m, string k, double d );
	vector< bool >   getVB( parameterMap m, string k, const vector< bool >&   d );
	vector< int >    getVI( parameterMap m, string k, const vector< int >&    d );
	vector< double > getV(  parameterMap m, string k, const vector< double >& d );
	vector< string > getVS( parameterMap m, string k, const vector< string >& d );
	vector< vector< double > > getM( parameterMap m, string k, const vector< vector< double > >& d );
	
	// Getter functions with no defaults
	bool   getB( parameterMap m, string k );
	int    getI( parameterMap m, string k );
	string getS( parameterMap m, string k );
	double getD( parameterMap m, string k );
	vector< bool >   getVB( parameterMap m, string k );
	vector< int >    getVI( parameterMap m, string k );
	vector< double > getV(  parameterMap m, string k );
	vector< string > getVS( parameterMap m, string k );
	vector< vector< double > > getM( parameterMap m, string k );
	parameterMap getMap( parameterMap m, string k );
	
	// Setters 
	void set( parameterMap& m, string k, bool   v );
	void set( parameterMap& m, string k, int    v );
	void set( parameterMap& m, string k, double v );
	void set( parameterMap& m, string k, string v );
	void set( parameterMap& m, string k, char*  v );
	void set( parameterMap& m, string k, vector< bool >   v );
	void set( parameterMap& m, string k, vector< int >    v );
	void set( parameterMap& m, string k, vector< double > v );
	void set( parameterMap& m, string k, vector< string > v );
	void set( parameterMap& m, string k, vector< vector< double > > v );
	void set( parameterMap& m, string k, parameterMap v );
	
	//-------------------------------------------------
	// Read parameterMap from a file
	//-------------------------------------------------
	void ReadParsFromFile(parameterMap& apars, string inFileName);
	void ReadParsFromFile(parameterMap& apars, char *inFileName);
	
	//-------------------------------------------------
	// Write parameterMap to a file
	//-------------------------------------------------
	void WriteParsToFile(parameterMap& apars, string outFileName);
	
	void PrintPars(parameterMap&);
	
		//! Stream insertion operator for parameterMaps
	//ostream& operator<<( ostream& out_stream, const parameterMap& p );
	
	//! Stream extraction operator for parameterMaps
	//istream& operator>>( istream& in_stream, parameterMap& p );

	
};

//! Stream insertion operator for parameterMaps
ostream& operator<<( ostream& out_stream, const parameterMap& p );

//! Stream extraction operator for parameterMaps
istream& operator>>( istream& in_stream, parameterMap& p );

#endif
