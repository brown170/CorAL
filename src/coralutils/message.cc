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
//---------------------------------------------------
/*! \class CMessage
\brief A nice interface to the message system

This class is an interface to the CMessage class. It provides an easy to work with interface similar to cerr, and adds a few other features.
*/
//-----------------------------------------------------

#include "message.h"
#include <iostream>
#include <stdio.h>

/*!
This makes the message class a singleton
*/
CMessage * CMessage::getCMessage() 
{
  static CMessage __me;
  return & __me;
}

/*!
This is the default constructer and when it is called the output is sent to the standard place, cerr.
*/
CMessage::CMessage(): show_info(false), show_warn(true), _constructorArgument("") {}

/*! This constructor will cause the messages to be sent with the options.*/
CMessage::CMessage(string msg): show_info(false), show_warn(true), _constructorArgument(msg) {}

/*! These are the operator overloads that allow you to use this class like cerr*/
CMessage& CMessage::operator<<(string str)
{
  _queue = _queue + str;
  return *this;
}

CMessage& CMessage::operator<<(const char* str)
{
  _queue = _queue + str;
  return *this;
}

CMessage& CMessage::operator<<(int in)
{
  char tmp[100];
  if(_hex) {
    sprintf(tmp,"%x",in);
  }else{
    sprintf(tmp,"%i",in);
  }
  _queue = _queue + tmp;
  return *this;
}

CMessage& CMessage::operator<<(long in)
{
  char tmp[100];
  if(_hex) {
    sprintf(tmp,"%lx",in);
  }else{
    sprintf(tmp,"%li",in);
  }
  _queue = _queue + tmp;
  return *this;
}

CMessage& CMessage::operator<<(double in)
{
  char tmp[100];
  sprintf(tmp,"%f",in);
  _queue = _queue + tmp;
  return *this;
}

CMessage& CMessage::operator<<(level lev)
{
  _level = lev;
  return *this;
}

CMessage& CMessage::operator<<(cntl test)
{
  switch(test)
    {
    case hex:
      _hex=true;
      break;
    case dec:
      _hex=false;
      break;
    case endm:
      if(_constructorArgument==""){
        if(_level==info && CMessage::show_info)    {cout<<"Info::"<<_queue<<endl;}
        if(_level==warning && CMessage::show_warn) {cerr<<"Warning::"<<_queue<<endl;}
        if(_level==severe)  {cerr<<"SevereError::"<<_queue<<endl;}
        if(_level==fatal)   {cerr<<"FatalError::"<<_queue<<endl;}
      }else{
        cout<<"other strings not implemented yet"<<endl;
      }
      _queue="";
      break;
    default:
      cerr<<"ERROR: Problem in message passing system "<<__FILE__<<":"<<__LINE__<<endl;
    }
  return *this;
}

