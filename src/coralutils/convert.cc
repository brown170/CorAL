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
#include "convert.h"
#include "message.h"


//!Removes the leading blanks in a string.
string remove_leading_blanks( const string& instring )
{
    if(( instring.size() <= 0 ) || (instring.find(' ')==string::npos) ) return( instring );
    int nspace=instring.find_first_not_of(' ');
    return instring.substr( nspace, instring.size()-nspace+1);
}

//!Removes the trailing blanks in a string.
string remove_trailing_blanks( const string& instring )
{
    if(( instring.size() <= 0 ) || (instring.find(' ')==string::npos) ) return( instring );
    int nspace=instring.find_last_not_of(' ');
    return instring.substr( 0, nspace+1);
}

//!Removes all extra blanks from a string.
string remove_extra_blanks( const string& instring )
{ 
    vector<string> split_string = split(instring);
    vector<string> no_blanks;
    for (vector<string>::iterator itr=split_string.begin(); itr!=split_string.end(); ++itr){
        if( (*itr).size() != 0 ) no_blanks.push_back(*itr);
    }
    return join(no_blanks," ");
}

//! Removes all blanks from a string.
string remove_all_blanks( const string& instring )
{
    vector<string> split_string = split(instring);
    vector<string> no_blanks;
    for (vector<string>::iterator itr=split_string.begin(); itr!=split_string.end(); ++itr){
        if( (*itr).size() != 0 ) no_blanks.push_back(*itr);
    }
    return join(no_blanks,"");
}

//! A function to remove everything past a comment
string remove_comments( const string& instring, const string& comment_string ){
    if (instring.find(comment_string)==string::npos) return instring;
    int new_string_length = instring.find(comment_string);
    return instring.substr( 0, new_string_length );
}

//! splits a string at all occurances of "pattern" and returns 
//! the substrings in a list 
vector< string > split( const string& s, const string pattern )
{
    vector< string > result;
    int begin_slice = 0;
    int end_slice = 0;
    int slice_length = 0;
    while ( s.find( pattern, begin_slice ) != string::npos )
    {
        end_slice = s.find( pattern, begin_slice )+1;
        slice_length = end_slice - begin_slice - 1;
        if ( slice_length > 0 )
            result.push_back( s.substr( begin_slice, slice_length ) );
        begin_slice = end_slice+pattern.size()-1;
    }
    result.push_back( s.substr( begin_slice ) );
    return result;
}

//! joins a list of strings using "pattern" as glue
string join( const vector< string >& ls, const string& pattern )
{
    string result="";
    for ( vector<string>::const_iterator itr = ls.begin(); itr != ls.end(); ++itr )
    {
        vector<string>::const_iterator next_itr = itr;
        ++next_itr;
        result += *itr;
        if ( next_itr != ls.end() ) result += pattern;
    }
    return result;
}

//! A function that lower-cases strings
string tolower( const string& s )
{
    string t = s;
    int length = t.length();
    for (int i = 0; i < length; i++) t[i] = tolower(t[i]);
    return t;
}

//! i copies of a string x pasted together
string operator*( int i, const string& x ){
    string xTimes("");
    for (int j=0;j<i;++j) xTimes+=x;
    return xTimes;
}

//! center a string in a field of a certain length
string center( const string& stuff, int len ){
    int numSpaces = len-stuff.size();
    numSpaces/=2;
    string aSpace(" ");
    string output(numSpaces*aSpace+stuff+numSpaces*aSpace);
    if (output.size()==static_cast<unsigned int>(len)) return output; else return output+" ";
}

//! left justify a string in a field of a certain length
string ljust( const string& stuff, int len ){
    int numSpaces = len-stuff.size();
    string aSpace(" ");
    return stuff+numSpaces*aSpace;
}

//! right justify a string in a field of a certain length
string rjust( const string& stuff, int len ){
    int numSpaces = len-stuff.size();
    string aSpace(" ");
    return numSpaces*aSpace+stuff;
}


//! A function that removes all instances of a string
string remove( const string& instring, const string& target ){
    return replace( instring, target, "");
}

//! A function that replaces all instances of one string with another
string replace( const string& instring, const string& old_pattern, const string& new_pattern ){
    return join( split( instring, old_pattern ), new_pattern );
}

//! A function to slice a string, using Python-like indexing
string slice( const string& instring, int istart, int ifinish ){
    int ibegin, iend;
    if (istart<0) { ibegin = instring.size()+istart;}  
    if (ifinish<0) { iend = instring.size()+ifinish;}
    return instring.substr( ibegin, iend-ibegin );
}


