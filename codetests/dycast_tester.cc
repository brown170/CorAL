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
#include "../include/locobject.h"
#include <string>
#include <iostream>
#include <vector>

class A: public CObject{
    public:
        A(void){}
        virtual ~A(void){}
        string method(void){return "method from A";}
};

class B: public A{
    public:
        string methodB(void){return "new method from B";}
};

class C: public B{
};

class D{
    public:
        string methodD(void){return "method from extra class D";}
};

class E: public C, public D{
    public:
        string method(void){return "override method from A";}
        string methodB(void){return "override method from B";}
};


int main(void){
    A aObj;
    B bObj;
    E eObj;
    cout << "\nObject Tests:"<<endl<<"============="<<endl;
    cout << "aObj.method(): "<<aObj.method()<<endl;    
    cout << "bObj.method(): "<<bObj.method()<<endl;    
    cout << "bObj.methodB(): "<<bObj.methodB()<<endl;    
    cout << "eObj.method(): "<<eObj.method()<<endl;    
    cout << "eObj.methodB(): "<<eObj.methodB()<<endl;    
    cout << "eObj.methodD(): "<<eObj.methodD()<<endl;    

    A* ePtr(&eObj);
    cout << "\nPointer Tests:"<<endl<<"=============="<<endl;
    cout << "ePtr has type A*, but points to an E"<<endl;
    cout << "ePtr->method(): "<< ePtr->method()<< endl;
    cout << "ePtr->methodB(): doesn't work, A's don't have methodB" << endl; 
    cout << "dynamic_cast<B*>(ePtr)->methodB(): "<<dynamic_cast<B*>(ePtr)->methodB()<<endl;
    cout << "dynamic_cast<E*>(ePtr)->method(): "<<dynamic_cast<E*>(ePtr)->method()<<endl;
    cout << "dynamic_cast<E*>(ePtr)->methodD(): "<<dynamic_cast<E*>(ePtr)->methodD()<<endl;

    CLocator<A> loE(new E);
    cout << "\nLocator Tests:"<<endl<<"=============="<<endl;
    cout << "loE has type CLocator<A>, but points to an E"<<endl;
    cout << "loE->method(): "<< loE->method()<< endl;
    cout << "loE->methodB(): doesn't work, A's don't have methodB" << endl; 
    cout << "dynamic_cast<B*>(&(*loE))->methodB(): "<<dynamic_cast<B*>(&(*loE))->methodB()<<endl;
    cout << "dynamic_cast<E*>(&(*loE))->method(): "<<dynamic_cast<E*>(&(*loE))->method()<<endl;
    cout << "dynamic_cast<E*>(&(*loE))->methodD(): "<<dynamic_cast<E*>(&(*loE))->methodD()<<endl<<endl;

    vector< CLocator<A> > vecOfE(0);
    cout << "\nVector of Locator Tests:"<<endl<<"========================"<<endl;
    for (unsigned int i=0;i<3;++i){
        CLocator<A> tmp(new E);
        vecOfE.push_back(tmp);
        cout << "vecOfE["<<i<<"]->method(): "<< vecOfE[i]->method()<< endl;
        cout << "vecOfE["<<i<<"]->methodB(): doesn't work, A's don't have methodB" << endl; 
        cout << "dynamic_cast<B*>(&(*(vecOfE["<<i<<"])))->methodB(): "<<dynamic_cast<B*>(&(*(vecOfE[i])))->methodB()<<endl;
        cout << "dynamic_cast<E*>(&(*(vecOfE["<<i<<"])))->method(): "<<dynamic_cast<E*>(&(*(vecOfE[i])))->method()<<endl;
        cout << "dynamic_cast<E*>(&(*(vecOfE["<<i<<"])))->methodD(): "<<dynamic_cast<E*>(&(*(vecOfE[i])))->methodD()<<endl;
        cout << endl;
    }
    
    return true;
}
