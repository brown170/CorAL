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
#ifndef CONVERT_UTILITIES
#define CONVERT_UTILITIES

#include <string>
#include <iostream>
#include <vector>

using namespace std;

// prototypes of conversion functions

//! A function that removes the leading blanks of a string
string remove_leading_blanks( const string& instring );

//! A function that removes the trailing blanks of a string
string remove_trailing_blanks( const string& instring );

//! A function that removes extra blanks in a string
string remove_extra_blanks( const string& instring );

//! A function that removes all blanks in a string
string remove_all_blanks( const string& instring );

//! A function that removes all instances of a string
string remove( const string& instring, const string& target );

//! A function that replaces all instances of one string with another
string replace( const string& instring, const string& old_pattern, const string& new_pattern );

//! A function to slice a string
string slice( const string& instring, int istart, int ifinish );

//! A function to remove everything past a comment
string remove_comments( const string& instring, const string& comment_string );

//! splits a string at all occurances of "pattern" and returns 
//! the substrings in a list 
vector< string > split( const string& s, const string pattern=" " );

//! joins a list of strings using "pattern" as glue
string join( const vector< string >& ls, const string& pattern="," );

//! A function that lower-cases strings
string tolower( const string& s );

//! i copies of a string x pasted together
string operator*( int i, const string& x );

//! center a string in a field of a certain length
string center( const string& stuff, int len );

//! left justify a string in a field of a certain length
string ljust( const string& stuff, int len );

//! right justify a string in a field of a certain length
string rjust( const string& stuff, int len );

#endif

