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
#ifndef NEW_YLMEXPANSION_H
#define NEW_YLMEXPANSION_H

#include "objects3d.h"
#include "objects1d.h"
#include "sf.h"
#include "cheezyparser.h"
#include "utils.h"
#include "constants.h"
#include <complex>
#include <map>
#include <sstream>
#include <fstream>

typedef complex<double> double_complex;

#define EXTENSION ".coral"

using namespace std;

//----------------------------------------------
///  Simple classes to term indices & basis functions
/// Base class, defines l
//----------------------------------------------
class CHarmonicBasisFunction{
public:
    int l;
    CHarmonicBasisFunction(int _l=0): l(_l){}
    CHarmonicBasisFunction(const CHarmonicBasisFunction& other):l(other.l){}
    virtual ~CHarmonicBasisFunction(void){}
    bool operator==(const CHarmonicBasisFunction& other) const{return other.l==l;}
    bool operator<(const CHarmonicBasisFunction& other) const{return l<other.l;}
    virtual double getValue(double theta, double phi) const=0;
    virtual double getValue(double x, double y, double z) const=0;
    string termName(void) const{string s;stringstream ss; ss<<"_"<<l; ss>>s; return s;}
};

//----------------------------------------------
/// For Ylm's have l, m and real/imag flag
//----------------------------------------------
class CSphericalHarmonicBasisFunction: public CHarmonicBasisFunction{
public:
    int l;
    int m;
    bool realpart;
    CSphericalHarmonicBasisFunction(int _l=0, int _m=0, bool _re=true): 
        l(_l), m(_m), realpart(_re){}
    CSphericalHarmonicBasisFunction(const CSphericalHarmonicBasisFunction& other):
        l(other.l),m(other.m),realpart(other.realpart){}
    virtual ~CSphericalHarmonicBasisFunction(void){} 
    bool operator==(const CSphericalHarmonicBasisFunction& other) const
        {return (other.l==l)&&(other.m==m)&&(other.realpart==realpart);}
    bool operator<(const CSphericalHarmonicBasisFunction& other) const{
        if (l<other.l) return true; 
        else if (m<other.m) return true; 
        else if (realpart==other.realpart) return false; 
        else return other.realpart;
    }
    double getValue(double theta, double phi) const
        {if (realpart) return SpherHarmonics::ReYlm(l,m,theta,phi); else return SpherHarmonics::ImYlm(l,m,theta,phi);}
    double getValue(double x, double y, double z) const
        {if (realpart) return SpherHarmonics::ReYlm(l,m,x,y,z); else return SpherHarmonics::ImYlm(l,m,x,y,z);}
    string termName(void) const
        {string s;stringstream ss; ss<<"_"<<l<<"_"<<m<<"_"; if (realpart) ss<<"re"; else ss<<"im"; ss>>s; return s;}
};

//----------------------------------------------
/// Cartesian harmonics have lx, ly and lz 
//----------------------------------------------
class CCartesianHarmonicBasisFunction: public CHarmonicBasisFunction{
public:
    int lx, ly, lz;
    CCHCalc cartHarm;
    CCartesianHarmonicBasisFunction(int _lx=0, int _ly=0, int _lz=0): 
        lx(_lx),ly(_ly),lz(_lz){}
    CCartesianHarmonicBasisFunction(const CCartesianHarmonicBasisFunction& other):
        lx(other.lx),ly(other.ly),lz(other.lz){}
    virtual ~CCartesianHarmonicBasisFunction(void){} 
    bool operator==(const CCartesianHarmonicBasisFunction& other) const
        {return (other.lx==lx)&&(other.ly==ly)&&(other.lz==lz);}
    bool operator<(const CCartesianHarmonicBasisFunction& other) const{
        int sum1 = lx+ly+lz, sum2=other.lx+other.ly+other.lz;
        if (sum1!=sum2) return (sum1<sum2); 
        if (lx!=other.lx) return (lx<other.lx); 
        if (ly!=other.ly) return (ly<other.ly); 
        else return (lz<other.lz); 
    }
    double getValue(double theta, double phi) const{
        double nx = sin(theta), ny = nx*sin(phi), nz = cos(theta);
        nx*=cos(phi);
        return const_cast<CCartesianHarmonicBasisFunction*>(this)->cartHarm.GetAFromE(lx,ly,lz,nx,ny,nz);
    }
    double getValue(double x, double y, double z) const
        {   
            double r = sqrt( x*x+y*y+z*z );
            if ( r < 1e-10 ) return const_cast<CCartesianHarmonicBasisFunction*>(this)->cartHarm.GetAFromE(lx,ly,lz,0.,0.,0.);
            return const_cast<CCartesianHarmonicBasisFunction*>(this)->cartHarm.GetAFromE(lx,ly,lz,x/r,y/r,z/r);
        }
    string termName(void) const
        {string s;stringstream ss; ss<<"_"<<lx<<"_"<<ly<<"_"<<lz; ss>>s; return s;}
};

//----------------------------------------------
///  Simple interface for all objects that are 
///  expansions, either Spherical or Cartesian 
///  It is templated, so everything is here in header
///
///  Don't use it directly! Inherit from it.  See
///  the CSphericalHarmonicExpansion class below for example.
///
///  TBas Should inherit from CHarmonicBasisFunction and TObj Should inherit from CObject1d
//----------------------------------------------
template< class TBas, class TObj >  
class CHarmonicExpansion : public CObject3d, public map< TBas, TObj > {

public:

    string storage_directory;
    string file_prefix;
    int lmax;
    bool skip_odd_l;
    typedef typename CHarmonicExpansion< TBas, TObj >::iterator term_iterator;
    typedef typename CHarmonicExpansion< TBas, TObj >::const_iterator const_term_iterator;

    // Constructors, Destructors
    CHarmonicExpansion(void): storage_directory("."), file_prefix("term"), lmax(0), skip_odd_l(false) {}
    CHarmonicExpansion(int _lmax, bool _skip, string _sd, string _fp = "term" ): 
        storage_directory(_sd), file_prefix(_fp), lmax(_lmax), skip_odd_l(_skip) {}
    virtual ~CHarmonicExpansion(void){}
        
    // ----------- read -----------------
    //! Read from parameterMap
    bool Read(const parameterMap& s){
        bool result(CObject3d::Read(s));
        lmax = parameter::getI(s,"lmax",0);
        if (s.find("skip_odd_l")==s.end())
            MESSAGE << "Looking for skip_odd_l flag, make sure you have reset this in any class that inherits from CHarmonicExpansion"<<ENDM_WARN;
        skip_odd_l = parameter::getB(s,"skip_odd_l",true);
        if ((skip_odd_l)&&(IS_ODD(lmax))) {lmax=lmax-1;}
        // Read in the filename & storage directory
        storage_directory = parameter::getS(s,"storage_directory",storage_directory);
        file_prefix = parameter::getS(s,"file_prefix",file_prefix);
        return result;
    }

    // ----------- write -----------------
    //! Write to parameterMap 
    bool Write(parameterMap& s){
        bool result=CObject3d::Write(s);
        parameter::set(s,"lmax",lmax);
        parameter::set(s,"skip_odd_l",skip_odd_l);
        parameter::set(s,"file_prefix",file_prefix);
        parameter::set(s,"storage_directory",storage_directory);
        return result;
    }
    
    // ----------- writeTerms -----------------
    //! Write terms to disk
    virtual bool writeTerms( void ){
        bool result=true;
        for ( term_iterator it=CHarmonicExpansion< TBas, TObj >::begin(); it!=CHarmonicExpansion< TBas, TObj >::end(); ++it ){
            if (keepTerm(it)) {
                parameterMap outParams;
                it->second.Write(outParams);
                string filename = storage_directory+"/"+file_prefix+(it->first.termName())+EXTENSION;
                ofstream oFile( filename.c_str() );
                if (!oFile) throw MESSAGE << "Writing term "<<it->first.termName()<<" to "<< filename <<" failed"<<ENDM_SEVERE;
                result = result&&(oFile<<outParams);
            }
        }
        return result;
    }

    // ----------- fillTerms ---------------
    //! Fill out the list with pre-initialized, but empty, terms
    virtual bool fillTerms( void )=0;

    // ----------- keepTerm -----------------
    //! test if we should keep this term
    bool keepTerm( const_term_iterator it ) const
        {return (it->first.l<=lmax) && !( IS_ODD(it->first.l)&&(skip_odd_l) );}

    // ----------- getValueSphr -----------------
    //! get value at specific location, in spherical coordinates
    double getValueSphr(double r, double theta, double phi) const{
        double ans=0.0;
        for ( const_term_iterator it=CHarmonicExpansion< TBas, TObj >::begin(); it!=CHarmonicExpansion< TBas, TObj >::end(); ++it ){
            if (keepTerm(it)) 
                ans += SQRTFOURPI*(it->second.getValue(r))*(it->first.getValue(r,theta,phi));
        }
        return ans;
    }

    // ----------- getErrorSphr -----------------
    //!  get uncertainty on value at specific location, in spherical coordinates
    double getErrorSphr(double r, double theta, double phi) const{
        double ans=0.0;
        for ( const_term_iterator it=CHarmonicExpansion< TBas, TObj >::begin(); it!=CHarmonicExpansion< TBas, TObj >::end(); ++it ){
            if (keepTerm(it)) 
                ans+=SQRTFOURPI*(it->second.getError(r))*(it->first.getValue(r,theta,phi));
        }
        return ans;
    }


    // ----------- getValueCart -----------------
    //!  get value at specific location, in Cartesian coordinates
    double getValueCart(double rS, double rO, double rL) const{
        double r=sqrt(rS*rS+rO*rO+rL*rL),ans=0.0;
        for ( const_term_iterator it=CHarmonicExpansion< TBas, TObj >::begin(); it!=CHarmonicExpansion< TBas, TObj >::end(); ++it ){
            if (keepTerm(it)) 
                ans+=SQRTFOURPI*(it->second.getValue(r))*(it->first.getValue(rS,rO,rL));
        }
        return ans;
    }

    // ----------- getErrorCart -----------------
    //!  get uncertainty on value at specific location, in Cartesian coordinates
    double getErrorCart(double rS, double rO, double rL) const{
        double r=sqrt(rS*rS+rO*rO+rL*rL),ans=0.0;
        for ( const_term_iterator it=CHarmonicExpansion< TBas, TObj >::begin(); it!=CHarmonicExpansion< TBas, TObj >::end(); ++it){
            if (keepTerm(it)) 
                ans+=SQRTFOURPI*(it->second.getError(r))*(it->first.getValue(rS,rO,rL));
        }
        return ans;
    }

};


//----------------------------------------------
///  For Spherical Harmonic expansion, add some
///  simple functions to interface 
//----------------------------------------------
template< class TObj >
class CSphericalHarmonicExpansion: public CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >{

public:

    typedef typename CSphericalHarmonicExpansion< TObj >::iterator term_iterator;
    typedef typename CSphericalHarmonicExpansion< TObj >::const_iterator const_term_iterator;

    // Constructors, Destructors
    CSphericalHarmonicExpansion(void){}
    CSphericalHarmonicExpansion(int _lmax, bool _skip, string _sd=".", string _fp="term"): 
         CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >(_lmax,_skip,_sd,_fp){}
    ~CSphericalHarmonicExpansion(void){}

    // ----------- readTerms -----------------
    //! Read terms from disk
    virtual bool readTerms( void ){
        bool result = false;
        int lstep=1;  if (this->skip_odd_l) lstep=2;
        for (int l=0; l<=this->lmax; l+=lstep){
            for (int m=0; m<=l; ++m){
                bool reim=false;
                for (int im=0; im<2; ++im){
                    reim=!reim;
                    if (!(m==0 && !reim)) {
                        CSphericalHarmonicBasisFunction key(l,m,reim);
                        cout << "    Reading term "<<key.termName()<<flush;
                        string filename = this->storage_directory+"/"+this->file_prefix+(key.termName())+EXTENSION;
                        ifstream iFile( filename.c_str() );
                        parameterMap inParams;
                        cout <<" .."<<flush;
                        if (iFile) { 
                            iFile >> inParams;
                            result = true;
                            TObj obj;
                            cout <<".."<<flush;
                            obj.Read(inParams);
                            cout <<".."<<flush;
                            this->insert(make_pair(key,obj));
                            cout <<" OK"<<endl;
                        } else {
                            MESSAGE << "Could not load term w/ filename "<<filename<<ENDM_WARN;
                        }
                    }
                }
            }
        }
        return result;
    }
    
    // ----------- fillTerms ---------------
    //! Fill out the list with pre-initialized, but empty, terms
    bool fillTerms( void ){
        bool result = false;
        int lstep=1;  if (this->skip_odd_l) lstep=2;
        for (int l=0; l<=this->lmax; l+=lstep){
            for (int m=0; m<=l; ++m){
                bool reim=false;
                for (int im=0; im<2; ++im){
                    reim=!reim;
                    if (!(m==0 && !reim)) {
                        CSphericalHarmonicBasisFunction key(l,m,reim);
                        TObj obj;
                        obj.l=l;
                        obj.m=m;
                        obj.realpart=reim;
                        this->insert(make_pair(key,obj));
                    }
                }
            }
        }
        return result;
    }


    // ----------- getValueSphr -----------------
    //! get complex value of term w/ lm, specific r
    double_complex getValueSphr(int l, int m, double r){
        const_term_iterator riter=find(l,abs(m),true);
        const_term_iterator iiter=find(l,abs(m),false);
        double repart(0.), impart(0.);
        if (riter!=CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >::end()) repart=riter->second.getValue(r);
        if (iiter!=CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >::end()) impart=iiter->second.getValue(r);
        if (m>=0) {
            return double_complex(repart, impart);
        } else {
            return static_cast<double_complex>(NEGONE_TO_THE(abs(m)))*double_complex(repart, -impart);
        }
    }

    // ----------- getErrorSphr -----------------
    //! get complex uncertainty of term w/ lm, specific r
    double_complex getErrorSphr(int l, int m, double r){  
        const_term_iterator riter=find(l,abs(m),true);
        const_term_iterator iiter=find(l,abs(m),false);
        double repart(0.), impart(0.);
        if (riter!=CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >::end()) repart=riter->second.getError(r);
        if (iiter!=CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >::end()) impart=iiter->second.getError(r);
        if (m>=0) {
            return double_complex(repart, impart);
        } else {
            return static_cast<double_complex>(NEGONE_TO_THE(abs(m)))*double_complex(repart, -impart);
        }
    }

    // ----------- find -----------------
    //! find a term in the map, return an iterator
    term_iterator find(int l, int m, bool rpart){
        return CHarmonicExpansion< CSphericalHarmonicBasisFunction, TObj >::find(CSphericalHarmonicBasisFunction(l,m,rpart));
    }
    
    // ------------ getItem ----------------
    //! must pick favorite interface.  Note:
    //! this is going to throw an error
    //! if it can't find the term
    TObj& getItem(int l, int m, bool rpart)
        {return (*this)[CSphericalHarmonicBasisFunction(l,m,rpart)];}
    TObj& operator()(int l, int m, bool rpart)
        {return (*this)[CSphericalHarmonicBasisFunction(l,m,rpart)];}
    TObj getItem(int l, int m, bool rpart) const
        {return (*this)[CSphericalHarmonicBasisFunction(l,m,rpart)];}
    TObj operator()(int l, int m, bool rpart) const
        {return (*this)[CSphericalHarmonicBasisFunction(l,m,rpart)];}
};

#endif
