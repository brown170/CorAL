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

#if !defined LOCOBJECT_H
#define      LOCOBJECT_H

#define __USING_RTTI__

#include "parametermap.h"
//#include "fileparser.h"
#include "message.h"
#include <string>
#include <map>

//-------------------------------------------------------------------------
// Declarations for class CObject.
//-------------------------------------------------------------------------

class CObject {
    friend class CLocatorBase;
private:
    long instcount;
    static long objcount;
    string objtype;
protected:
    CObject(void) : instcount(0L), objtype("") {}
    CObject(const CObject& m) : instcount(0L), objtype("") {}
    virtual ~CObject(void) {}
public:
    virtual CObject& CopyState(const CObject& m);

    long InstCount(void) const {return instcount;}
    static long ObjCount(void) {return objcount;}

    virtual CObject* Clone(void) const;

    virtual void setObjType(string t){objtype = t;}
    virtual string getObjType(void){return objtype;}
    
    // fancy pointer detailed report
    void ReportObjectPars(ostream& strm) const;

    // read/write to CComandOptions map (must override these)
    virtual bool Read(const parameterMap& s){return false;}
    virtual bool Write(parameterMap& s){return false;}

//    // read/write to CCommand (DO NOT OVERRIDE THESE)
//    virtual bool Read(const CCommand& c);
//    virtual bool Write(CCommand& c);

};

//-------------------------------------------------------------------------
// Declarations for class CLocatorBase.
//-------------------------------------------------------------------------

class CLocatorBase {
protected:
    CObject* pObj;
public:
    CLocatorBase(void) : pObj(NULL) {}
    explicit CLocatorBase(CObject* _pObj) : pObj(NULL)
        {AttachObject(_pObj);}
    CLocatorBase(const CLocatorBase& lo);
    ~CLocatorBase(void) {Dispose();}
    CLocatorBase& operator=(const CLocatorBase& lo);
    CLocatorBase& AttachObject(CObject* _pObj);
    CLocatorBase& AttachClonedObject(
        const CLocatorBase& lo, bool copystate);
    CObject& operator*(void)
        {if (pObj == NULL) ReportError();  return *pObj;}
    const CObject& operator*(void) const
        {if (pObj == NULL) ReportError();  return *pObj;}
    CObject* operator->(void)
        {if (pObj == NULL) ReportError();  return pObj;}
    const CObject* operator->(void) const
        {if (pObj == NULL) ReportError();  return pObj;}
    long InstanceCount(void) const {return pObj ? pObj->instcount : 0L;}
    bool IsValid(void) const {return pObj != NULL;}
    void DetachObject(void) {AttachObject(NULL);}
    void Report(ostream& strm) const;
private:
    void Dispose(void);
protected:
    void ReportError(void) const;
};

inline ostream& operator<<(ostream& strm, const CLocatorBase& lo)
    {lo.Report(strm);return strm;}

//-------------------------------------------------------------------------
// Declarations for class template CLocator.
//-------------------------------------------------------------------------

template <class T>
class CLocator : public CLocatorBase {
public:
    CLocator(void) : CLocatorBase() {}
    explicit CLocator(CObject* _pObj) : CLocatorBase(_pObj) {}
    CLocator(const CLocator<T>& lo) : CLocatorBase(lo) {}
    CLocator<T>& operator=(const CLocator<T>& lo)
        {
            *static_cast<CLocatorBase*>(this) = static_cast<const CLocatorBase&>(lo);
            return *this;
        }
    CLocator<T>& AttachObject(T* _pObj)
        {return static_cast<CLocator<T>&>(CLocatorBase::AttachObject(_pObj));}
    CLocator<T>& AttachClonedObject(const CLocator<T>& lo, bool copystate)
        {return static_cast<CLocator<T>&>(CLocatorBase::AttachClonedObject(lo, copystate));}
    operator T*(void)
        {if (pObj == NULL) ReportError();return static_cast<T*>(pObj);}
    operator const T*(void) const
        {if (pObj == NULL) ReportError();return static_cast<const T*>(pObj);}
    T& operator*(void)
        {if (pObj == NULL) ReportError();return *static_cast<T*>(pObj);}
    const T& operator*(void) const
        {if (pObj == NULL) ReportError();return *static_cast<const T*>(pObj);}
    T* operator->(void)
        {if (pObj == NULL) ReportError();return static_cast<T*>(pObj);}
    const T* operator->(void) const
        {if (pObj == NULL) ReportError();return static_cast<const T*>(pObj);}
    T* operator()(void)
        {if (pObj == NULL) ReportError();return static_cast<T*>(pObj);}
    const T* operator()(void) const
        {if (pObj == NULL) ReportError();return static_cast<const T*>(pObj);}
};

template <class T>
bool operator< (const CLocator<T>& l, const CLocator<T>& r)
{
    const T* plObj = l;
    const T* prObj = r;

    if (plObj == NULL) return true;
    if (prObj == NULL) return false;

    return *plObj < *prObj;
}

template <class T>
bool operator==(const CLocator<T>& l, const CLocator<T>& r)
{
    const T* plObj = l;
    const T* prObj = r;

    if ((plObj == NULL) && (prObj == NULL)) return true;
    if ((plObj == NULL) || (prObj == NULL)) return false;

    return *plObj == *prObj;
}

//-------------------------------------------------------------------------
// Declarations for class CObjectRegistry.
//-------------------------------------------------------------------------

class CObjectRegistry : public map<string, CLocatorBase> {
public:
    void Register(const CLocatorBase& lo, const string& ID);
    bool Unregister(const string& ID);
    CObject* Create(const string& ID) const;
    void ClearRegistry(void) {clear();}
    void Report(ostream& strm) const;
};

inline ostream& operator<<(ostream& strm, const CObjectRegistry& m)
{
    m.Report(strm);
    return strm;
}

//-------------------------------------------------------------------------
// Macros to register CObject objects in CObjectRegistry reg
// and to retrieve them attached to Locator objects.
// ID is an identifier string, enclosed in double quotes.
//-------------------------------------------------------------------------

#define REGISTER_OBJECT(reg, CObjDerived, ID) \
    { \
        string id = ID; \
        CLocator<CObjDerived> lo(new CObjDerived); \
        reg.Register(lo, id); \
        lo->setObjType(id); \
    }

#define CREATE_OBJECT(reg, CObjTarg, lo, ID) \
    CLocator<CObjTarg> lo; \
    { \
        string id = ID; \
        CObject* pObj = reg.Create(id); \
        PTR_CAST(CObjTarg*, pObjTarg, pObj); \
        lo.AttachObject(pObjTarg); \
        lo->setObjType(id); \
    }

//-------------------------------------------------------------------------
// Macros to register CObject objects in GlobalObjectRegistry
// and to retrieve them attached to Locator objects.
// The identifier is the class name.
//-------------------------------------------------------------------------

#if !defined LOCOBJECT_CPP
    extern CObjectRegistry GlobalObjectRegistry;
#endif

#define REGISTER_GLOBAL_OBJECT(CObj) \
    { \
        CLocator<CObj> lo(new CObj); \
        GlobalObjectRegistry.Register(lo, #CObj); \
    }

#define UNREGISTER_GLOBAL_OBJECT(CObj) \
    GlobalObjectRegistry.Unregister(#CObj);

#define CREATE_GLOBAL_OBJECT(CObjTarg, lo, CObj) \
    CLocator<CObjTarg> lo; \
    { \
        CObject* pObj = GlobalObjectRegistry.Create(#CObj); \
        PTR_CAST(CObjTarg*, pObjTarg, pObj); \
        lo.AttachObject(pObjTarg); \
    }

#endif
//-------------------------------------------------------------------------
// Macros to allow dynamic casting of CLocators
//-------------------------------------------------------------------------
#define DYNAMIC_CAST(CObj,ClassName) dynamic_cast<ClassName*>(&(*CObj))

//-------------------------------------------------------------------------
// Macros for downcasting polymorphic types.
//-------------------------------------------------------------------------

#if defined __USING_RTTI__

	#define PTR_CAST(T, DestPtr, SrcPtr) \
		T DestPtr = dynamic_cast<T>(SrcPtr); \
		if (DestPtr == NULL) {\
			MESSAGE<<"Dynamic cast of "#SrcPtr" to type "#T" failed."<<ENDM_FATAL; \
            exit(-1);\
		} 
        
	#define REF_CAST(T, DestRef, SrcRef) \
		try { \
			T DestRef = dynamic_cast<T>(SrcRef); \
		} \
		catch (...) { \
			MESSAGE<<"Dynamic cast of "#SrcRef" to type "#T" failed."<<ENDM_FATAL; \
            exit(-1);\
		} \
		T DestRef = static_cast<T>(SrcRef);

#else

	#define PTR_CAST(T, DestPtr, SrcPtr) \
		T DestPtr = static_cast<T>(SrcPtr);

	#define REF_CAST(T, DestRef, SrcRef) \
		T DestRef = static_cast<T>(SrcRef);

#endif

//-------------------------------------------------------------------------
// Macros to allow casting away const-ness of member data from stuff 
// pointed to by CLocators
//-------------------------------------------------------------------------
#define CONST_CAST(Data,DataType) *const_cast<DataType*>(&Data)
