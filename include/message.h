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
#ifndef _MESSAGE_
#define _MESSAGE_

#include <string>

using namespace std;

#define __ENDM_BASE__ CMessage::dec<<" [file:"<<__FILE__<<" line:"<<__LINE__<<" function:"<<__FUNCTION__<<"]"<<CMessage::endm

#define ENDM_INFO   CMessage::info<<__ENDM_BASE__
#define ENDM_WARN   CMessage::warning<<__ENDM_BASE__
#define ENDM_SEVERE CMessage::severe<<__ENDM_BASE__
#define ENDM_FATAL  CMessage::fatal<<__ENDM_BASE__

#define MESSAGE *CMessage::getCMessage()

class CMessage
{

 public:
 
  static CMessage * getCMessage(void);
  bool show_info;
  bool show_warn;

  enum level {info,warning,severe,fatal};
  enum cntl {hex,dec,endm};

  CMessage& operator<< (string);
  CMessage& operator<< (const char*);
  CMessage& operator<< (int);
  CMessage& operator<< (unsigned int tmp){*this<<(int)tmp; return *this;};
  CMessage& operator<< (long);
  CMessage& operator<< (double);

  CMessage& operator<< (cntl);
  CMessage& operator<< (level);

 private:  
 
  CMessage();
  CMessage(string);
  
  string _queue;
  level _level;
  string _constructorArgument;
  bool _hex;

};

#endif
