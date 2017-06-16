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
#ifndef EXPANDER_H
#define EXPANDER_H

#include <vector>
#include <complex>
#include <iostream>
#include "parametermap.h"
#include "message.h"
#include "histogram1d.h"
#include "histogram3d.h"
#include "harmonic_expansion.h"
#include "sf.h"
#include "constants.h"
#include "utils.h"
#include "tnt_array1d.h"
#include "corr3d_histo.h"
#include "corr3d_ylm.h"
#include "sou3d_histo.h"
#include "sou3d_ylm.h"
#include "sou1d_histo.h"
#include "linalg.h"
#include "integratevec.h"

using namespace TNT;
using namespace parameter;
using namespace std;

#define CONST_CAST(Data,DataType) *const_cast<DataType*>(&Data)

// Main converter functions, use these!  
// They set the member data not present in the CHistogram3d and
// CObjectYlmExpansion base classes.  Although the C[..]Ftn3dSphr
// classes don't specify the nature of the 1d terms, these functions
// create/require only C[..]Ftn1dHisto terms.  These functions do 
// test for them too. 
CCorrFtn3dSphr    expand( const CCorrFtn3dHisto& corrin, const parameterMap m );
CCorrFtn3dHisto   recombine( const CCorrFtn3dSphr& corrin, const parameterMap m  );
CSourceFtn3dSphr<CSourceFtn1dHisto>  expand( const CSourceFtn3dHisto& souin, const parameterMap m  );
CSourceFtn3dHisto recombine( const CSourceFtn3dSphr<CSourceFtn1dHisto>& souin, const parameterMap m  );

template< class THisto1d >
class CCart2SphrExpander {
    
public:
    
    //-------------------- Constructor ------------------------
    // Constructor
    CCart2SphrExpander( void ): 
    num_int_vec( 9 ), l_list( num_int_vec ), 
    m_list( num_int_vec ), reim_list( num_int_vec ), 
    mtx_list( num_int_vec ){} 
    
    //-------------------- Member data... ------------------------
    // For easy access, these are public, so try not to break them
    int lmax;
    int ndata;
    int nr;
    bool generate_full_covmtx;
    // Flags to tell whether data is "flip-symmetric"
    bool xflipflag;
    bool yflipflag;
    bool zflipflag;
    double symfactor;
    
    //-------------------- Read ------------------------
    //! Read from parameter map
    bool Read( const parameterMap& m ){
        CONST_CAST(lmax,int) = parameter::getI(m,"lmax",0);
        CONST_CAST(num_int_vec,int) = parameter::getI(m,"num_int_vec",9);
        CONST_CAST(generate_full_covmtx,bool) = parameter::getB(m,"generate_full_covmtx",false);
        return true;
    }
    
    //-------------------- Write ------------------------
    //! Write to parameter map
    bool Write( parameterMap& m ){
        parameter::set(m,"lmax",lmax);
        parameter::set(m,"num_int_vec",num_int_vec);
        parameter::set(m,"generate_full_covmtx",generate_full_covmtx);
        return true;
    }
    
    // ---------------- convert ------------------
    //! Convert to Spherical binning
    bool convert( const CHistogram3d& sou, CSphericalHarmonicExpansion< THisto1d >& targ, bool skip_odd_l=true ){
        // Make tmp object, with correct memory allocation.  
        // At the end or this routine, we'll copy its guts into targ.
        // We have to do it this way to ensure that the correct type of terms
        // are copied into the final targ.  
        targ.lmax=lmax;
        targ.skip_odd_l=skip_odd_l;
        // Initialize targ and the CSphericalHarmonicExpander
        double dr=(sou.dx+sou.dy+sou.dz)/3.0;
        double rmax=std::min(std::min(sou.rightBinEdgeX(sou.nx),sou.rightBinEdgeY(sou.ny)),sou.rightBinEdgeZ(sou.nz));
        nr = iround<double>((rmax-dr/2.0)/dr);
        int lstep=1;  if (skip_odd_l) lstep=2;
        for (int l=0; l<=lmax; l+=lstep){
            for (int m=0; m<=l; ++m){
                { // for all l,m there is a real part
                    CSphericalHarmonicBasisFunction key(l,m,true);
                    THisto1d obj;
                    obj.l=l;
                    obj.m=m;
                    obj.realpart=true;
                    obj.redim(nr);
                    obj.setFixedWidthBins(dr,dr/2);
                    targ.insert(make_pair(key,obj));
                }
                if (m!=0) { // for m!=0, there is an imaginary part too
                    CSphericalHarmonicBasisFunction key(l,m,false);
                    THisto1d obj;
                    obj.l=l;
                    obj.m=m;
                    obj.realpart=false;
                    obj.redim(nr);
                    obj.setFixedWidthBins(dr,dr/2);
                    targ.insert(make_pair(key,obj));
                }
            }
        }
        // Create and set up converter matrices
        ndata = sou.ndata;
        setSymmetryFactorAndFlags( sou, skip_odd_l );
        cout << "    Some CCart2SphrExpander properties"<<endl;
        cout << "        dr: " << dr<< "; rmax: " <<rmax << endl;
        cout << "        ndata: " << ndata<< "; size: " <<targ.size()<< "; size/num_int_vec+1: " << targ.size()/num_int_vec+1<< endl;
        cout << "        num_int_vec = " << num_int_vec<<"; nr = "<<nr<<endl;
        cout << endl;
        // Compute terms, in groups of num_int_vec
        typename CSphericalHarmonicExpansion< THisto1d>::iterator chunk_start=targ.begin();
        typename CSphericalHarmonicExpansion< THisto1d>::iterator element;
        int ielement;
        int ichunk=0;
        while ( chunk_start!=targ.end() ){
            resizeLists(0);
            ielement=0;
            for ( element = chunk_start; (element!=targ.end())&&(ielement<num_int_vec); ++element ) {
                l_list.push_back(element->first.l);
                m_list.push_back(element->first.m);
                reim_list.push_back(element->first.realpart);
                ++ielement;
            } 
            int chunk_size=l_list.size();
            // process integration list
            MESSAGE << "   Building matrix group #" << ichunk << "  chunk_size: "<<chunk_size<<ENDM_INFO; 
            setUpMtxList( chunk_size, nr, ndata );
            buildCart2SphrMtx( sou, targ.begin()->second );
            MESSAGE << "   Building  CSphericalHarmonicExpansion ..." << ENDM_INFO;
            MESSAGE << "     ielement  l  m  Real" << ENDM_INFO;
            ielement=0;
            for ( element=chunk_start; (element!=targ.end())&&(ielement<num_int_vec); ++element ) {
                MESSAGE << "     " << ielement << "         " 
                        << l_list[ielement] <<"  " << m_list[ielement] <<"  " 
                        << reim_list[ielement] ;
                element->second.data = mtx_list[ielement]*sou.data;
                MESSAGE << "       data done ... " ;
                element->second.covmtx = ErrorsToCovmtx(mtx_list[ielement],sou.uncert);
                MESSAGE << "covmtx done ... " ;
                element->second.syncUncert();
                MESSAGE << "uncert done ... " << ENDM_INFO;
                CONST_CAST(element->second.covmtx_is_active,bool)=true;
                ++ielement;
            }
            chunk_start=element;
            ++ichunk;
        }
        return true;
    }
        
    // ---------------- convert ------------------
    //! Convert back to Cartesian binning
    bool convert( const CSphericalHarmonicExpansion< THisto1d >& sou, CHistogram3d& targ ){    
        nr=sou.begin()->second.ndata; // all the terms have same radial bins, but we'll check
        double dr=sou.begin()->second.binWidth(0);
        double rmax=sou.begin()->second.leftBinEdge(nr);
        int Nx=nr;
        double x0=dr/2.0;
        if (!sou.skip_odd_l) {Nx=2*nr;x0=-rmax+dr/2.;}
        // Check all the terms in the sou, to make sure they are all compatible with this 
        // function.  We require the terms to all be CHistogram1d's and to have the same
        // fixed width binning.
        for (typename CSphericalHarmonicExpansion<THisto1d>::const_iterator it=sou.begin();it!=sou.end();++it){
            if (!it->second.fixed_width_bins){
                MESSAGE << "CSphericalHarmonicExpansion term "
                        << it->first.l<<" "<<it->first.m<<" "<<it->first.realpart
                        << " does not have fixed width bins" << ENDM_FATAL;
                return false;
            }
            if (it->second.ndata!=nr){
                MESSAGE << "CSphericalHarmonicExpansion term "
                        << it->first.l<<" "<<it->first.m<<" "<<it->first.realpart
                        <<" does not have same number of bins as other terms: " 
                        << nr << " vs. "<<it->second.ndata<<ENDM_FATAL;
                return false;
            }
            if (abs(it->second.binWidth(0)-dr)>0.001*dr){
                MESSAGE << "CSphericalHarmonicExpansion term "
                        << it->first.l<<" "<<it->first.m<<" "<<it->first.realpart
                        <<" does not have same bin size as other terms: " 
                        << it->second.binWidth(0) << " vs. "<<dr <<ENDM_FATAL;
                return false;
            }
        }
        // Make tmp object, with correct memory allocation.  
        // and properly initialized
        targ.nx=Nx;
        targ.ny=2*nr;
        targ.nz=2*nr;
        targ.dx=dr;
        targ.dy=dr;
        targ.dz=dr;
        targ.xoffset=x0;
        targ.yoffset=-rmax+dr/2.;
        targ.zoffset=-rmax+dr/2.;
        targ.redim(targ.nx*targ.ny*targ.nz);
        targ.ixzero=targ.findBin(targ.dx/2.,targ.dx,targ.xoffset);
        targ.iyzero=targ.findBin(targ.dy/2.,targ.dy,targ.yoffset);
        targ.izzero=targ.findBin(targ.dz/2.,targ.dz,targ.zoffset);
        // Initialize the CSphericalHarmonicExpander internal data
        lmax=sou.lmax;
        ndata=targ.ndata;
        // Create and set up converter matrices
        setSymmetryFactorAndFlags( targ, sou.skip_odd_l );
        cout << "   Unexpanding ... " << endl;
        // setup temp arrays
        Array1D<double> data(ndata,0.0);
        Array1D<double> error(ndata,0.0);
        Array1D<bool> dud(ndata,false);
        // Compute terms, in chunks w/ chunk_size <= num_int_vec
        typename CSphericalHarmonicExpansion<THisto1d>::const_iterator chunk_start=sou.begin();
        typename CSphericalHarmonicExpansion<THisto1d>::const_iterator element;
        int ielement;
        int ichunk=0;
        while ( chunk_start!=sou.end() ){
            // build integration list of size num_int_vec
            resizeLists(0);
            ielement=0;
            for ( element=chunk_start; (element!=sou.end())&&(ielement<num_int_vec); ++element){
                l_list.push_back(element->first.l);
                m_list.push_back(element->first.m);
                reim_list.push_back(element->first.realpart);
                ++ielement;
            }             
            int chunk_size=l_list.size();
            // process integration list
            MESSAGE << "   Building matrix group #" << ichunk << "  chunk_size: "<<chunk_size<<ENDM_INFO; 
            setUpMtxList( chunk_size, ndata, nr );
            buildSphr2CartMtx( sou.begin()->second, targ );
            MESSAGE << "   Building CHistogram3d ..." << ENDM_INFO;
            MESSAGE << "     ielement  l  m  Real" << ENDM_INFO;
            ielement=0;
            for ( element=chunk_start; (element!=sou.end())&&(ielement<num_int_vec); ++element ) {
                MESSAGE << "     " << ielement << "            " 
                        << l_list[ielement] <<"  " << m_list[ielement] <<"  " 
                        << reim_list[ielement]  <<"    "<<ENDM_INFO;
                // unexpand each term
                Array1D<double> tdat, terr;
                tdat=mtx_list[ielement]*(element->second.data);
                terr=CovmtxToErrors(mtx_list[ielement],element->second.covmtx);
                // must add in the unexpanded terms
                for (int j=0;j<ndata;++j){
                    data[j] += tdat[j];
                    if (generate_full_covmtx) {
                        // nothing yet
                    } else {
                        error[j] = sqrt(error[j]*error[j]+terr[j]*terr[j]);
                    }
                }
                ++ielement;
            }
            chunk_start=element;
            ++ichunk;
        }        
        // setup the dudlist
        for (int i=0;i<targ.nx;++i){
            for (int j=0;j<targ.ny;++j){
                for (int k=0;k<targ.nz;++k){
                    // any bin that is not completely covered by the radial grid is considered a dud
                    double len1,len2,len3,len4,len5,len6,len7,len8, maxlen;
                    len1=length( targ.leftBinEdgeX(i)/2.0,  targ.leftBinEdgeY(j)/2.0,  targ.leftBinEdgeZ(k)/2.0  );
                    len2=length( targ.rightBinEdgeX(i)/2.0, targ.leftBinEdgeY(j)/2.0,  targ.leftBinEdgeZ(k)/2.0  );
                    len3=length( targ.leftBinEdgeX(i)/2.0,  targ.rightBinEdgeY(j)/2.0, targ.leftBinEdgeZ(k)/2.0  );
                    len4=length( targ.leftBinEdgeX(i)/2.0,  targ.leftBinEdgeY(j)/2.0,  targ.rightBinEdgeZ(k)/2.0 );
                    len5=length( targ.rightBinEdgeX(i)/2.0, targ.rightBinEdgeY(j)/2.0, targ.leftBinEdgeZ(k)/2.0  );
                    len6=length( targ.rightBinEdgeX(i)/2.0, targ.leftBinEdgeY(j)/2.0,  targ.rightBinEdgeZ(k)/2.0 );
                    len7=length( targ.leftBinEdgeX(i)/2.0,  targ.rightBinEdgeY(j)/2.0, targ.rightBinEdgeZ(k)/2.0 );
                    len8=length( targ.rightBinEdgeX(i)/2.0, targ.rightBinEdgeY(j)/2.0, targ.rightBinEdgeZ(k)/2.0 );
                    maxlen=MAXF(MAXF(MAXF(len1,len2),MAXF(len3,len4)),MAXF(MAXF(len5,len6),MAXF(len7,len8)));
                    if (maxlen>sou.begin()->second.rightBinEdge(nr)) dud[targ.whatIndex(i,j,k)]=true;
                }
            }
        }
        // OK, for all the points that are clearly duds, let's reset those correlation points
        for (int i=0;i<ndata;++i){if (dud[i]) {data[i]=1.0;error[i]=1.0;}}        
        // set the final results
        targ.data = data;
        targ.uncert = error;
        targ.covmtx_is_active=false;
        return true;
    } 
        
private:
        
    // ----------------- resizeLists ----------------
    // Pre-run setup functions
    void resizeLists( int nSize ){
        l_list.resize( nSize ); 
        m_list.resize( nSize );
        reim_list.resize( nSize );
        mtx_list.resize( nSize );
    }
    
    // ----------------- setSymmetryFactorAndFlags ----------------
    // Pre-run setup functions
    void setSymmetryFactorAndFlags( const CHistogram3d& sou, bool skip_odd_l ){
        symfactor = 1.0;
        xflipflag = false;
        yflipflag = false;
        zflipflag = false;
        if (skip_odd_l) {
            // if any of following true, then no q<0 in data in that dimension
            if (sou.ixzero==0) {xflipflag=true;symfactor*=2.0;} 
            if (sou.iyzero==0) {yflipflag=true;symfactor*=2.0;}
            if (sou.izzero==0) {zflipflag=true;symfactor*=2.0;}
        } 
        MESSAGE << "   Zero Bins are (" << sou.ixzero << ", " 
            << sou.iyzero << ", " << sou.izzero
            << ") so setting CExpander symmetry factor to "
            <<symfactor<<ENDM_INFO;
    }
    
    // ----------------- setUpMtxList ----------------
    // Simple Matrix Ops.
    void setUpMtxList( int nSize, int nrows, int ncols ){
        mtx_list.resize(0);
        for (int i=0;i<nSize;++i){
            Array2D<double> temp(nrows,ncols,0.0);
            mtx_list.push_back(temp);
        }
    }
    
    // ------------------ buildSphr2CartMtx -------------------
    void buildSphr2CartMtx( const CHistogram1d& souTerm, const CHistogram3d& targ ){
        double factor=SQRTFOURPI/targ.binVolume();
        double qmin, qmax;
        int l,m;
        bool reim;
        
        // integrator for sphr->cart matrix
        CIntegrateVector J;
        J.SetNDim(3);
        J.SetNumFunc(mtx_list.size());
        J.SetMaxPts(10000);
        
        cout << "     Radial bin counter: "<<flush;
        for (int n=0;n<nr;++n){
            rlo=souTerm.leftBinEdge(n);
            rhi=souTerm.rightBinEdge(n);
            
            cout << n <<" "<<flush;
            
            int iSmin=targ.ixzero;
            int iSmax=targ.nx;
            
            for (int iS=iSmin;iS<iSmax;++iS){
                
                // get edges of this bin
                xlo=targ.midBinX(iS)-targ.dx/2.0;
                xhi=targ.midBinX(iS)+targ.dx/2.0;
                
                // figure out limits of next loop given this bin    
                int iOmin=targ.iyzero;
                int iOmax=targ.ny;
                
                for (int iO=iOmin;iO<iOmax;++iO){
                    
                    // get edges of this bin
                    ylo=targ.midBinY(iO)-targ.dy/2.0;
                    yhi=targ.midBinY(iO)+targ.dy/2.0;
                    
                    // figure out limits of next loop given this bin    
                    int iLmin=targ.izzero;
                    int iLmax=targ.nz;
                    
                    for (int iL=iLmin;iL<iLmax;++iL){
                        
                        // get edges of this bin
                        zlo=targ.midBinZ(iL)-targ.dz/2.0;
                        zhi=targ.midBinZ(iL)+targ.dz/2.0;
                        
                        qmax=length(xhi,yhi,zhi);
                        qmin=length(xlo,ylo,zlo);
                        
                        if (((qmin<=souTerm.midBin(n))&&(qmax>=souTerm.midBin(n)))||
                            ((qmin<=rhi)&&(qmax>=rhi))||
                            ((qmin<=rlo)&&(qmax>=rlo))){
                            
                            // fixed limits
                            J.SetLimits(0,xlo,xhi);
                            J.SetLimits(1,ylo,yhi);
                            J.SetLimits(2,zlo,zhi);
                            
                            // Integrate!
                            J.Compute(this,sphr2CartMtxIntegrand);// for sphr->cart
                            
                            // Unpack into mtxlist
                            for (unsigned int ilist=0;ilist<mtx_list.size();++ilist){
                                l=l_list[ilist];
                                m=m_list[ilist];
                                reim=reim_list[ilist];
                                
                                // precompute flipped indices
                                int a    =targ.whatIndex( iS,                iO,                iL                );
                                int a_y  =targ.whatIndex( iS,                targ.flipBinY(iO), iL                );
                                int a_z  =targ.whatIndex( iS,                iO,                targ.flipBinZ(iL) );
                                int a_yz =targ.whatIndex( iS,                targ.flipBinY(iO), targ.flipBinZ(iL) );
                                // these may not be used if x < 0 flipped onto x > 0
                                int a_x  =targ.whatIndex( targ.flipBinX(iS), iO,                iL                );
                                int a_xy =targ.whatIndex( targ.flipBinX(iS), targ.flipBinY(iO), iL                );
                                int a_xz =targ.whatIndex( targ.flipBinX(iS), iO,                targ.flipBinZ(iL) );
                                int a_xyz=targ.whatIndex( targ.flipBinX(iS), targ.flipBinY(iO), targ.flipBinZ(iL) );
                                
                                double mtxelement = factor*J.GetResults(ilist);
                                if (reim) { // REAL PART....
                                            // use symmetry to set the elements as needed
                                    mtx_list[ilist][a][n]     =                      mtxelement;
                                    if (!yflipflag)                           mtx_list[ilist][a_y][n]   =                      mtxelement;
                                    if (!zflipflag)                           mtx_list[ilist][a_z][n]   = NEGONE_TO_THE(l+m) * mtxelement;
                                    if (!(yflipflag&&zflipflag))              mtx_list[ilist][a_yz][n]  = NEGONE_TO_THE(l+m) * mtxelement;
                                    if (!xflipflag)                           mtx_list[ilist][a_x][n]   = NEGONE_TO_THE(m)   * mtxelement;
                                    if (!(xflipflag&&yflipflag))              mtx_list[ilist][a_xy][n]  = NEGONE_TO_THE(m)   * mtxelement;
                                    if (!(xflipflag&&zflipflag))              mtx_list[ilist][a_xz][n]  = NEGONE_TO_THE(l)   * mtxelement;
                                    if (!(xflipflag&&(yflipflag&&zflipflag))) mtx_list[ilist][a_xyz][n] = NEGONE_TO_THE(l)   * mtxelement;
                                    
                                } else { // IMAGINARY PART....
                                         // use symmetry to set the elements as needed
                                    mtx_list[ilist][a][n]     =                       mtxelement;
                                    if (!yflipflag)                           mtx_list[ilist][a_y][n]   = -1.0                * mtxelement;
                                    if (!zflipflag)                           mtx_list[ilist][a_z][n]   =  NEGONE_TO_THE(l+m) * mtxelement;
                                    if (!(yflipflag&&zflipflag))              mtx_list[ilist][a_yz][n]  = -NEGONE_TO_THE(l+m) * mtxelement;
                                    if (!xflipflag)                           mtx_list[ilist][a_x][n]   = -NEGONE_TO_THE(m)   * mtxelement;
                                    if (!(xflipflag&&yflipflag))              mtx_list[ilist][a_xy][n]  =  NEGONE_TO_THE(m)   * mtxelement;
                                    if (!(xflipflag&&zflipflag))              mtx_list[ilist][a_xz][n]  = -NEGONE_TO_THE(l)   * mtxelement;
                                    if (!(xflipflag&&(yflipflag&&zflipflag))) mtx_list[ilist][a_xyz][n] =  NEGONE_TO_THE(l)   * mtxelement;
                                }
                            }
                        }
                    }
                }
            }
        }
        cout <<endl;
    } 
    
    
    // ------------------ buildCart2SphrMtx -------------------
    void buildCart2SphrMtx( const CHistogram3d& sou, const CHistogram1d& targTerm ){
        double factor=1.0/SQRTFOURPI/targTerm.binWidth(0); //Should have fixed width bins
        double qmin, qmax;
        int l,m;
        bool reim;
        
        // integrator for cart->sphr matrix
        CIntegrateVector J;
        J.SetNDim(3);
        J.SetNumFunc(mtx_list.size());
        J.SetMaxPts(10000);
        
        cout << "     Radial bin counter: "<<flush;
        for (int n=0;n<nr;n++){
            rlo=targTerm.leftBinEdge(n);
            rhi=targTerm.rightBinEdge(n);
            
            cout << n <<" "<<flush;
            
            int iSmin=sou.ixzero;
            int iSmax=sou.nx;
            
            for (int iS=iSmin;iS<iSmax;iS++){
                
                // get edges of this bin
                xlo=sou.midBinX(iS)-sou.dx/2.0;
                xhi=sou.midBinX(iS)+sou.dx/2.0;
                
                // figure out limits of next loop given this bin    
                int iOmin=sou.iyzero;
                int iOmax=sou.ny;
                
                for (int iO=iOmin;iO<iOmax;iO++){
                    
                    // get edges of this bin
                    ylo=sou.midBinY(iO)-sou.dy/2.0;
                    yhi=sou.midBinY(iO)+sou.dy/2.0;
                    
                    // figure out limits of next loop given this bin    
                    int iLmin=sou.izzero;
                    int iLmax=sou.nz;
                    
                    for (int iL=iLmin;iL<iLmax;iL++){
                        
                        // get edges of this bin
                        zlo=sou.midBinZ(iL)-sou.dz/2.0;
                        zhi=sou.midBinZ(iL)+sou.dz/2.0;
                        
                        qmax=length(xhi,yhi,zhi);
                        qmin=length(xlo,ylo,zlo);
                        
                        if (((qmin<=targTerm.midBin(n))&&(qmax>=targTerm.midBin(n)))||
                            ((qmin<=rhi)&&(qmax>=rhi))||
                            ((qmin<=rlo)&&(qmax>=rlo))){
                            
                            // fixed limits
                            J.SetLimits(0,xlo,xhi);
                            J.SetLimits(1,ylo,yhi);
                            J.SetLimits(2,zlo,zhi);
                            
                            // Integrate!
                            J.Compute(this,cart2SphrMtxIntegrand);// for cart->sphr
                            
                            // Unpack into mtxlist
                            for (unsigned int ilist=0;ilist<mtx_list.size();ilist++){
                                l=l_list[ilist];
                                m=m_list[ilist];
                                reim=reim_list[ilist];
                                
                                // precompute flipped indices
                                int a    =sou.whatIndex( iS,               iO,               iL               );
                                int a_x  =sou.whatIndex( sou.flipBinX(iS), iO,               iL               );
                                int a_y  =sou.whatIndex( iS,               sou.flipBinY(iO), iL               );
                                int a_z  =sou.whatIndex( iS,               iO,               sou.flipBinZ(iL) );
                                int a_xy =sou.whatIndex( sou.flipBinX(iS), sou.flipBinY(iO), iL               );
                                int a_xz =sou.whatIndex( sou.flipBinX(iS), iO,               sou.flipBinZ(iL) );
                                int a_yz =sou.whatIndex( iS,               sou.flipBinY(iO), sou.flipBinZ(iL) );
                                int a_xyz=sou.whatIndex( sou.flipBinX(iS), sou.flipBinY(iO), sou.flipBinZ(iL) );
                                
                                // Because you only get a flipped index if it exists (otw., get back
                                // original index) and because the matrices are zeroed, we will
                                // add the matrix elements below. 
                                //  
                                // What happens is that, if there is no symmetry to worry about in a 
                                // particular direction, then we will need both the index and the 
                                // flipped index to fill out the matrix.  We'll just use Ylm symmetries 
                                // to generate those from the matrix in octant 1 and add them to zero.
                                // 
                                // On the other hand, if there is a symmetry in a particular direction, 
                                // when you flip, you get back the original index, so the adding operation
                                // will add the flipped matrix element into the unflipped one.
                                double mtxelement = factor*J.GetResults(ilist);
                                if (reim) { // REAL PART....
                                            // use symmetry to set the elements
                                    mtx_list[ilist][n][a]     +=                      mtxelement;
                                    mtx_list[ilist][n][a_x]   += NEGONE_TO_THE(m)   * mtxelement;
                                    mtx_list[ilist][n][a_y]   +=                      mtxelement;
                                    mtx_list[ilist][n][a_z]   += NEGONE_TO_THE(l+m) * mtxelement;
                                    mtx_list[ilist][n][a_xy]  += NEGONE_TO_THE(m)   * mtxelement;
                                    mtx_list[ilist][n][a_xz]  += NEGONE_TO_THE(l)   * mtxelement;
                                    mtx_list[ilist][n][a_yz]  += NEGONE_TO_THE(l+m) * mtxelement;
                                    mtx_list[ilist][n][a_xyz] += NEGONE_TO_THE(l)   * mtxelement;
                                    
                                } else { // IMAGINARY PART....
                                         // use symmetry to set the elements
                                    mtx_list[ilist][n][a]     +=                       mtxelement;
                                    mtx_list[ilist][n][a_x]   += -NEGONE_TO_THE(m)   * mtxelement;
                                    mtx_list[ilist][n][a_y]   += -1.0                * mtxelement;
                                    mtx_list[ilist][n][a_z]   +=  NEGONE_TO_THE(l+m) * mtxelement;
                                    mtx_list[ilist][n][a_xy]  +=  NEGONE_TO_THE(m)   * mtxelement;
                                    mtx_list[ilist][n][a_xz]  += -NEGONE_TO_THE(l)   * mtxelement;
                                    mtx_list[ilist][n][a_yz]  += -NEGONE_TO_THE(l+m) * mtxelement;
                                    mtx_list[ilist][n][a_xyz] +=  NEGONE_TO_THE(l)   * mtxelement;
                                }
                                
                            }
                                
                        }
                        
                    }
                    
                }
                
            }
            
        }
        
        cout <<endl;
        
    }
    
    
    
    // ---------------- rInner, rOuter -------------------
    //! For setting the integration limits
    double rInner( double x, double y ){return sqrt(MAXF(0.,rlo*rlo-x*x-y*y));}
    double rOuter( double x, double y ){return sqrt(MAXF(0.,rhi*rhi-x*x-y*y));}
    
    //! ------------------ cart2SphrMtxIntegrand -------------------
    // The integrand cart->sphr
    static void cart2SphrMtxIntegrand(void* classptr, int* ndim, double* q, int* numfunc, double* f){
        CCart2SphrExpander<THisto1d> *cls = (CCart2SphrExpander<THisto1d>*)classptr;
        double qrad=cls->length(q[0],q[1],q[2]);
        double qfact;
        if (qrad<1e-6) qfact=0.0; else qfact=1.0/qrad/qrad;
        double_complex ylm;
        for (unsigned int il=0; il<cls->mtx_list.size(); ++il){
            int l=cls->l_list[il];
            int m=cls->m_list[il];
            bool reim=cls->reim_list[il];
            if ((qrad<=cls->rhi)&&(qrad>=cls->rlo)) {
                ylm=SpherHarmonics::Ylm(l,m,q[0],q[1],q[2]);
                if (reim) {f[il]=real(ylm)*qfact;}
                else      {f[il]=-imag(ylm)*qfact;}
            } else {
                f[il]=0.0;
            }
        }
    }
    
    //! ------------------ sphr2CartMtxIntegrand -------------------
    // The integrand sphr->cart
    static void sphr2CartMtxIntegrand(void* classptr, int* ndim, double* q, int* numfunc, double* f){
        CCart2SphrExpander<THisto1d> *cls = (CCart2SphrExpander<THisto1d>*)classptr;
        double qrad=cls->length(q[0],q[1],q[2]);
        double_complex ylm;
        for (unsigned int il=0; il<cls->mtx_list.size(); ++il){
            int l=cls->l_list[il];
            int m=cls->m_list[il];
            bool reim=cls->reim_list[il];
            if ((qrad<=cls->rhi)&&(qrad>=cls->rlo)) {
                ylm=SpherHarmonics::Ylm(l,m,q[0],q[1],q[2]);
                if (m==0) {
                    if (reim) {f[il]=real(ylm);}
                    else {f[il]=0.0;}
                } else {
                    if (reim) {f[il]=2.0*real(ylm);}
                    else {f[il]=-2.0*imag(ylm);}
                }
            } else {
                f[il]=0.0;
            }
        }
    }
    
    // --------------- length, lengthsqr ------------------
    //! To get length of (x, y, z) vector
    double length( double qS, double qO, double qL ){return sqrt(qS*qS+qO*qO+qL*qL);}
    double lengthsqr( double qS, double qO, double qL ){return (qS*qS+qO*qO+qL*qL);}
    
    // --------------- private member data ------------------
    // Integration list arrays & their dimensions
    int num_int_vec;
    vector< int > l_list;
    vector< int > m_list;
    vector< bool > reim_list;
    vector< Array2D< double > > mtx_list;
    // Data to initialize result object
    double xlo, xhi, ylo, yhi, zlo, zhi, rlo, rhi;
    
};

#endif
