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
#include "cheezyparser.h"
#include "message.h"

//-------------------------------------------------
// utility functions for stream parser
//-------------------------------------------------
string get_prepped_line( istream& s, const string& comment_string ){
    string line;
    getline(s,line);
    return remove_trailing_blanks( remove_comments( line, comment_string ) );    
}

//-------------------------------------------------
// Data and functions for supported types
//-------------------------------------------------
// prototypes for some complicated classes
vector<int> __ivec__(1,1);
vector<float> __fvec__(1,1.0);
vector<double> __dvec__(1,1.0);
vector<bool> __bvec__(1,0);
vector<string> __svec__(1,"");
vector< vector< double > > __dmat__(1,__dvec__);
parameterMap __pmap__;

string get_type_string( const type_info& tid ){
    if      ( tid == typeid(int(1))     ) return   "int"          ;
    else if ( tid == typeid(float(1))   ) return   "float"        ;
    else if ( tid == typeid(double(1))  ) return   "double"       ;
    else if ( tid == typeid(bool(1))    ) return   "bool"         ;
    else if ( tid == typeid(string("")) ) return   "string"       ;
    else if ( tid == typeid("")         ) return   "string"       ;
    else if ( tid == typeid(__ivec__)   ) return   "vector_int"   ;
    else if ( tid == typeid(__fvec__)   ) return   "vector_float" ;
    else if ( tid == typeid(__dvec__)   ) return   "vector_double";
    else if ( tid == typeid(__bvec__)   ) return   "vector_bool"  ;
    else if ( tid == typeid(__svec__)   ) return   "free_text"    ;
    else if ( tid == typeid(__pmap__)   ) return   "parameter_map";
    else if ( tid == typeid(__pmap__)   ) return   "instance"     ;
    else if ( tid == typeid(__svec__)   ) return   "vector_string";
    else if ( tid == typeid(__dmat__)   ) return   "matrix_double";
    cerr  << "Type "<<tid.name()<<" is not supported!" << endl;
    return "";
}

boost::any stoany( const string& type, const string& instring ){
    if      ( type == "int" )           return boost::any(simple_lexical_cast<int>(instring));
    else if ( type == "float" )         return boost::any(simple_lexical_cast<float>(instring));
    else if ( type == "double" )        return boost::any(simple_lexical_cast<double>(instring));
    else if ( type == "bool" )          return boost::any(simple_lexical_cast<bool>(instring));
    else if ( type == "string" )        return boost::any( instring );
    else if ( type == "vector_int" )    return boost::any(simple_vector_lexical_cast< int >(instring));
    else if ( type == "vector_float" )  return boost::any(simple_vector_lexical_cast< float >(instring));
    else if ( type == "vector_double" ) return boost::any(simple_vector_lexical_cast< double >(instring));
    else if ( type == "vector_bool" )   return boost::any(simple_vector_lexical_cast< bool >(instring));
    else if ( type == "free_text" )     return boost::any( split(instring,"\n") );
    else if ( type == "vector_string" ) return boost::any(simple_vector_lexical_cast< string >(instring));
    else if ( type == "instance" )      return boost::any(simple_lexical_cast< parameterMap >(instring));
    else if ( type == "parameter_map" ) return boost::any(simple_lexical_cast< parameterMap >(instring));
    else if ( type == "matrix_double" ) {
    	vector< string > rows = split( instring, "\n" );
	vector< vector< double > > result;
	for ( vector< string >::iterator it = rows.begin(); it!=rows.end(); ++it ) 
		result.push_back( simple_vector_lexical_cast< double >(*it) );
    	return boost::any(result);
    }
    cerr  << "Making a boost::any out of type "+type+" is not supported!" << endl;
    return boost::any(0);
}

//-------------------------------------------------
//  stream based i/o
//-------------------------------------------------
//! Input for stl vectors of stl strings.
istream& operator>>( istream& i, vector< string > V ){
    string line;
    while (i.good()) {
        getline(i,line);
        if (line.find(END_BLOCK) == string::npos) 
            V.push_back( line );
        else return i;
    }
    return i;
}

//! Output for stl vectors of stl strings.
ostream& operator<<(ostream& o, vector< string > v){
    vector< string >::iterator it;
    for ( it=v.begin(); it!=v.end(); ++it ){
        o << (*it);
        ++it; if (it!=v.end()) o << endl; --it;
    } 
    return o;
}

//! Stream extraction operator for parameterMaps
istream& operator>>( istream& in_stream, parameterMap& p ){
    while ( in_stream ){
        // get a new line 
        string line = get_prepped_line(in_stream,COMMENT_STRING);
	if ( line == "" ) continue;

        // see if we're already done
        if ( line.find(END_BLOCK) != string::npos ) {return in_stream;}
        // if one-liner, it is easy
        if ( line.find(BEGIN_BLOCK) == string::npos ) {
             CParserObject parsed_line(line);
	     p.setItem(parsed_line.key, parsed_line.value);
        // OK, we apparently found a data block of some sort
        } else {
             CParserObject parsed_first_line(line);
	     //If the next item is a parameter_map, we fill it iteritivly
	     if(parsed_first_line.type == "parameter_map" || (parsed_first_line.type == "instance")){
	       parameterMap nestedMap;
	       in_stream>>nestedMap;
	       p.setItem(parsed_first_line.key, nestedMap);
	     } else {
               vector< string > block;
               bool found_end(false);
	       while ( in_stream.good() && !found_end ) {
                  line = get_prepped_line( in_stream, COMMENT_STRING );
                  found_end = ( line.find(END_BLOCK) != string::npos );
                  if ( (line!="") && !found_end ) block.push_back(line);
	       }
	       CParserObject parsed_block( parsed_first_line.type, parsed_first_line.key, join(block,"\n") );
	       p.setItem(parsed_block.key, parsed_block.value);
	     }
        }                                                
    }
    return in_stream;
}

//! Stream insertion operator for parameterMaps
ostream& operator<<( ostream& out_stream, const parameterMap& p ){    
    parameterMap::const_iterator itr;
    for (itr=p.begin();itr!=p.end();++itr){
        string type = get_type_string(itr->second.type());
        string key = itr->first;
        boost::any value = itr->second;
        out_stream << "\n"<< type << " " << key << " " << flush;
        // Little one liner types
        if      ( type == "int" )           out_stream << boost::any_cast<int>(value);
        else if ( type == "float" )         out_stream << boost::any_cast<float>(value);
        else if ( type == "double" )        out_stream << boost::any_cast<double>(value);
        else if ( type == "bool" )          out_stream << boolalpha << boost::any_cast<bool>(value);
        else if ( type == "string" )        out_stream << boost::any_cast<string>(value);
        // Big data blocks
        else {
            out_stream  << BEGIN_BLOCK << flush;      
            if      ( type == "vector_int" )    out_stream << "\n" << boost::any_cast< vector< int > >(value);
            else if ( type == "vector_float" )  out_stream << "\n" << boost::any_cast< vector< float > >(value);
            else if ( type == "vector_double" ) out_stream << "\n" << boost::any_cast< vector< double > >(value);
            else if ( type == "vector_bool" )   out_stream << "\n" << boolalpha << boost::any_cast< vector< bool > >(value);
            else if ( type == "free_text" )     out_stream << "\n" << boost::any_cast< vector< string > >(value);
            else if ( type == "vector_string" ) out_stream << "\n" << boost::any_cast< vector< string > >(value);
            else if ( type == "matrix_double" ) out_stream << "\n" << boost::any_cast< vector< vector< double > > >(value);
            else if ( type == "instance" )      out_stream << boost::any_cast< parameterMap >(value);
            else if ( type == "parameter_map" ) out_stream << boost::any_cast< parameterMap >(value);
            out_stream  << "\n" << END_BLOCK << flush;      
        }
    }
    return out_stream;
}

//-------------------------------------------------
//  CParserObject functions
//-------------------------------------------------
CParserObject::CParserObject( const string& instring ){
    vector< string > split_line = split(instring," ");
    type = split_line[0];
    key  = split_line[1];
    //We must skip past type occurence before searching for key in case type contains key
    string rest = instring.substr(instring.find(key, instring.find(type) + type.size())+key.size());
    if ( remove_extra_blanks(rest) == BEGIN_BLOCK ) {
        value = boost::any(0);
        return;
    }
    value = boost::any(stoany(type,remove_extra_blanks(rest)));
}

