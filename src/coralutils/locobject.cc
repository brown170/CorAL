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
// Version as of 10/20/02.

#define LOCOBJECT_CPP

#include "locobject.h"
#include "message.h"

#include <iomanip>

long CObject::objcount = 0L;

CObjectRegistry GlobalObjectRegistry;

//-------------------------------------------------------------------------
// Definitions for class CObject.
//-------------------------------------------------------------------------

CObject& CObject::CopyState(const CObject& m)
{
    MESSAGE<<"CopyState function appropriate to object is not available."<<ENDM_WARN;
    return *this = m;    // Use of m avoids possible warning
}

void CObject::ReportObjectPars(ostream& strm) const
{
    strm << "this               " << this              << endl
         << "instcount          " << instcount         << endl
         << "objcount           " << objcount          << endl;
}

CObject* CObject::Clone(void) const
{
    MESSAGE<<"Clone function appropriate to object is not available."<<ENDM_WARN;
    return NULL;
}

//bool CObject::Read(const CCommand& c){
//    if (c.arguments[0] == objtype) {return Read(c.parameters);}
//    else {return false;}
//}

//bool CObject::Write(CCommand& c){
//    c.arguments.push_back(objtype);
//    Write(c.parameters);
//    return true;
//}

//-------------------------------------------------------------------------
// Definitions for class CLocatorBase.
//-------------------------------------------------------------------------

CLocatorBase::CLocatorBase(const CLocatorBase& lo)
{
    pObj = lo.pObj;
    if (pObj) pObj->instcount++;
}

CLocatorBase& CLocatorBase::operator=(const CLocatorBase& lo)
{
    if (this != &lo) {
        Dispose();
        pObj = lo.pObj;
        if (pObj) pObj->instcount++;
    }
    return *this;
}

CLocatorBase& CLocatorBase::AttachObject(CObject* _pObj)
{
    Dispose();
    if (_pObj != NULL) {
        pObj = _pObj;
        pObj->instcount = 1L;
        CObject::objcount++;
    }else {
        pObj = NULL;
    }
    return *this;
}

CLocatorBase& CLocatorBase::AttachClonedObject(const CLocatorBase& lo, bool copystate)
{
    CObject* pNewObj = lo->Clone();
    AttachObject(pNewObj);
    if (copystate) pObj->CopyState(*lo);
    return *this;
}

void CLocatorBase::Report(ostream& strm) const
{
    strm << setw(12) << "pObj"
         << setw(12) << "instcount"
         << setw(12) << "objcount"
         << endl
         << setw(12) <<  pObj
         << setw(12) <<  (pObj ? pObj->instcount : 0L)
         << setw(12) <<  CObject::objcount
         << endl;
}

void CLocatorBase::Dispose(void)
{
    if (pObj != NULL) {
        pObj->instcount--;
        if (pObj->instcount == 0L) {
            CObject::objcount--;
            delete pObj;
            pObj = NULL;
        }
    }
}

void CLocatorBase::ReportError(void) const
{
    MESSAGE<<"Attempt to access an object that is not present."<<ENDM_WARN;
}

//-------------------------------------------------------------------------
// Definitions for class CObjectRegistry.
//-------------------------------------------------------------------------

void CObjectRegistry::Register(const CLocatorBase& lo, const string& ID)
{
    if (!insert(value_type(ID, lo)).second){MESSAGE<<"Duplicate identifier "<<ENDM_FATAL;exit(-1);}
}

bool CObjectRegistry::Unregister(const string& ID)
{
    iterator it = find(ID);
    if (it == end()) {
        return false;
    } else {
        erase(it);
        return true;
    }
}

CObject* CObjectRegistry::Create(const string& ID) const
{
    const_iterator it = find(ID);
    if (it == end()){MESSAGE<< "Identifier " + ID + " does not exist."<<ENDM_FATAL;exit(-1);}
    CObject* result = (*it).second->Clone();
    result->setObjType(ID);
    return result;
}

void CObjectRegistry::Report(ostream& strm) const
    {for (const_iterator it = begin(); it != end(); ++it)strm << "        " << (*it).first << endl;}
