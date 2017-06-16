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
#ifndef STREAMPARSER_H
#define STREAMPARSER_H

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <typeinfo>
#include <map>
#include "parametermap.h"
#include "any.h"
#include "convert.h"

#define COMMENT_STRING "#"
#define BEGIN_BLOCK "{"
#define END_BLOCK "}"

using namespace std;


//-------------------------------------------------
// utility functions for stream parser
//-------------------------------------------------
//! Readies a line for parsing
string get_prepped_line( istream& s, const string& comment_string=COMMENT_STRING );

//! Wrapper around stringstream to allow "casting" of strings
template< typename Target >
Target simple_lexical_cast( const string& argu ){
    stringstream interpreter;
    Target result;
    interpreter << argu << boolalpha;
    interpreter >> result;
    return result;
}

//! Wrapper around stringstream to allow "casting" of vectors of strings
template< typename Target >
vector< Target > simple_vector_lexical_cast( const string& argu ){
    stringstream interpreter;
    vector< Target > result;
    Target buff;
    interpreter << argu << boolalpha;
    while ( interpreter >> buff ) result.push_back( buff );
    return result;
}

//-------------------------------------------------
// Data and functions for supported types
//-------------------------------------------------
string get_type_string( const type_info& tid );

boost::any stoany( const string& type, const string& instring );

//-------------------------------------------------
//  stream based i/o
//-------------------------------------------------
//! Specific output for vectors of strings
istream& operator>>( istream& i, vector< string > V );

//! Specific input for vectors of strings
ostream& operator<<(ostream& o, vector< string > v);

//! Templated output for stl vectors
template< typename ValueType >
ostream& operator<<(ostream& o, vector< ValueType > v){
    typename vector< ValueType >::iterator it;
    for ( it=v.begin(); it!=v.end(); ++it ){
        o << (*it);
        ++it; if (it!=v.end()) o << endl; --it;
    } 
    return o;
}

//! Templated input for stl vectors.  Needs to know a little bit about the syntax.
template< typename ValueType >
istream& operator>>( istream& i, vector< ValueType > V ){
    string line;
    while (i.good()) {
        getline(i,line);
        if (line.find(END_BLOCK) == string::npos) {
            V.push_back( simple_lexical_cast< ValueType >(line) );
        }
        else return i;
    }
    return i;
}

//! Templated output for stl vector of vectors
template< typename ValueType >
ostream& operator<<(ostream& o, vector< vector< ValueType > > v){
    typename vector< vector< ValueType >  >::iterator it;
    typename vector< ValueType >::iterator rit;
    for ( it=v.begin(); it!=v.end(); ++it ){
        for ( rit=it->begin(); rit!=it->end(); ++rit ){
            o << (*rit) << " ";
        }
        ++it; if (it!=v.end()) o << endl; --it;
    } 
    return o;
}

//! Templated input for stl vector of vectors.  Needs to know a little bit about the syntax.
template< typename ValueType >
istream& operator>>( istream& i, vector< vector< ValueType > > V ){
    string line;
    V.clear();
    while (i.good()) {
        getline(i,line);
        if (line.find(END_BLOCK) == string::npos) {
            vector< string > srow = split(line);
            vector< double > row;
            for ( 
                vector< string >::iterator buff = srow.begin();
                buff != srow.end(); 
                ++buff
            ) row.push_back( simple_lexical_cast< double >(*buff) );
            V.push_back(row);
        }
        else return i;
    }
    return i;
}



//-------------------------------------------------
//  CParserObject 
//-------------------------------------------------
//! class to simplify parsing
class CParserObject{
  public:
    string      type;
    string      key;
    boost::any  value;
    CParserObject( const string& instring );
    CParserObject( const string& intype, const string& inkey, const string& inblock):
        type(intype), key(inkey), value( stoany( intype, inblock ) ){}
    CParserObject( const string& inkey, const boost::any& initem): 
        type(get_type_string(initem.type())), key(inkey), value(initem){}
};

#endif
