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
#include "message.h"
#include "parametermap.h"
#include "constants.h"
#include "utils.h"
#include "linalg.h"
#include "lsqrinvert.h"
#include "bspline_imager1d.h"

using namespace TNT;
using namespace std;

#define __USE_GSL__

//#define __IMAGER_EPSILON 1e-10
#define __TYPICAL_SOURCE_SCALE__ 1e-4 //in fm^-3


//---------------- Read ----------------------------
bool CBasisSplineImager1d::Read( const parameterMap& m ){

    CBasisFuncImager1d::Read( m );
    qscale                      = parameter::getD(m,"qscale",50.);
//    spline_degree               = parameter::getI(m,"spline_degree",3);
    knot_tolerance              = parameter::getD(m,"knot_tolerance",1e-2);
    knot_init_scheme            = parameter::getS(m,"knot_init_scheme","default");
    optimize_knots              = parameter::getB(m,"optimize_knots",false);
    if (optimize_knots) {
        optimize_knots_max_iter                 = parameter::getB(m,"optimize_knots_max_iter",optimize_knots_max_iter);
        optimize_knots_minimizer_func_weight    = parameter::getB(m,"optimize_knots_minimizer_func_weight",optimize_knots_minimizer_func_weight);
        optimize_knots_coellesce_knot_tol       = parameter::getB(m,"optimize_knots_coellesce_knot_tol",optimize_knots_coellesce_knot_tol);
        optimize_knots_simplex_min_size         = parameter::getB(m,"optimize_knots_simplex_min_size",optimize_knots_simplex_min_size);
    }
    vector< double > empty_vec(0);
    if (knot_init_scheme == "user_defined_collocation_points") {
        colloc_pts = parameter::getV(m,"user_collocation_points",empty_vec);
    }
    if (knot_init_scheme == "user_defined_knots") {
        user_knots = parameter::getV(m,"user_knots",empty_vec);
    }
    return true;
}

//--------------- Write -----------------------------
// not used much except for diagnostics, so write everything!
bool CBasisSplineImager1d::Write( parameterMap& m ){

    CBasisFuncImager1d::Write( m );
    parameter::set(m,"qscale",qscale);
//    parameter::set(m,"spline_degree",spline_degree);
    parameter::set(m,"knot_tolerance",knot_tolerance);
    parameter::set(m,"knot_init_scheme",knot_init_scheme);
    parameter::set(m,"optimize_knots",optimize_knots);
    parameter::set(m,"optimize_knots_max_iter",optimize_knots_max_iter);
    parameter::set(m,"optimize_knots_minimizer_func_weight",optimize_knots_minimizer_func_weight);
    parameter::set(m,"optimize_knots_coellesce_knot_tol",optimize_knots_coellesce_knot_tol);
    parameter::set(m,"optimize_knots_simplex_min_size",optimize_knots_simplex_min_size);
    parameter::set(m,"user_collocation_points",colloc_pts);
    parameter::set(m,"user_knots",user_knots);
    return true;
}
    
//---------------- convertCorrelationToSource ----------------------------
//! function that manages the imaging
bool CBasisSplineImager1d::convertCorrelationToSource( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m, const CKernel* _kernelPtr ){

    CBasisSpline1d* souPtr = dynamic_cast<CBasisSpline1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisSpline1d"<<ENDM_FATAL;
    bool result = CBasisFuncImager1d::convertCorrelationToSource( corrin, souout, m, _kernelPtr );
    if (optimize_knots && result) do_optimal_knots( *souPtr ); 
    return result;
}

//---------------- convertSourceToCorrelation ----------------------------
//! function that manages the unimaging
bool CBasisSplineImager1d::convertSourceToCorrelation( const CSourceFtnBase& souin, CCorrFtn1dHisto& corrout, const parameterMap& m, const CKernel* _kernelPtr ){

    return CBasisFuncImager1d::convertSourceToCorrelation( souin, corrout, m, _kernelPtr );
}

//--------------------- imageit ---------------------------
//! Code to image a 1d source, with equality constraints
double CBasisSplineImager1d::imageit(CBasisFunctionExpansion1d& souout){

    CBasisSpline1d* souPtr = dynamic_cast<CBasisSpline1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisSpline1d"<<ENDM_FATAL;

    // Image using the generic code or optimize knots
    double chi2 = CBasisFuncImager1d::imageit( *souPtr );

    // Compute final summed relative error of source (to emulate old "optimal knot" behavior)
    double relerr=0.0;
    for (int i=0;i<souout.ndata;++i){relerr+=abs(souout.data[i]/sqrt(souout.covmtx[i][i]));}
    
    // Compute the "figure of merit" that we'll minimize when doing the optimal knots
    return chi2 + optimize_knots_minimizer_func_weight * relerr;
}

//-------------------- do_optimal_knots ----------------------------
double CBasisSplineImager1d::do_optimal_knots(CBasisSpline1d& souout){
#ifdef __USE_GSL__

    // don't have to set souWork because we already ran imageit once
    // souWork.CopyState(souout);

    /* variable declarations */
    cout << "  Performing knot optimization"<<endl;
    size_t nVarKnots = ndata_sou-(spline_degree+1);
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
    gsl_multimin_fminimizer *s = NULL;
    gsl_vector *ss, *x;
    gsl_multimin_function minex_func;
    size_t iter = 1, i;  // we already ran imageit once, so set iter to 1
    int status;
    double size;

    /* Initial vertex size vector */
    ss = gsl_vector_alloc (nVarKnots);

    /* Set all step sizes to ?? */
    double dxKnots=(rmax)/(nVarKnots+1);
    gsl_vector_set_all (ss, 0.75*dxKnots);

    /* Starting point */
    x = gsl_vector_alloc (nVarKnots);
    for (unsigned int i=0;i<nVarKnots;++i) gsl_vector_set(x, i, souout.knots[i+spline_degree+1]);
    cout << "    Variable knots initialized"<<endl;

    /* Initialize method */
    minex_func.f = &static_imageit;
    minex_func.n = nVarKnots;
    minex_func.params = (void *)this;
    s = gsl_multimin_fminimizer_alloc (T, nVarKnots);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);
    cout << "    Simplex minimizer initialized"<<endl;

    /* Iterate */
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
      
        if (status) break;

        size = gsl_multimin_fminimizer_size (s);
        status = gsl_multimin_test_size (size, 1e-2);

        if (status == GSL_SUCCESS) printf ("converged to minimum at\n");

        printf ("%5d ", int(iter));
        for (i = 0; i < nVarKnots; i++) printf ("%12.5e ", gsl_vector_get (s->x, i));
        
        printf ("f() = %7.3f size = %.3f\n", s->fval, size);
    } while (status == GSL_CONTINUE && iter < 200);
    
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free (s);

    souout.CopyState(souWork);
    return (double)(status);

#else

    return imageit(souout);

#endif
}

//-------------------- static_imageit ----------------------------
//! Wrapper to imageit that we can use to do optimization
//! using GSL.  Also, recomputes losource's knots using
//! the variable knots in variableKnots, handed to this
//! routine by the gsl minimizer
double CBasisSplineImager1d::static_imageit(const gsl_vector *variableKnots, void* classPtr){

#ifdef __USE_GSL__
    CBasisSplineImager1d *cls = (CBasisSplineImager1d*)classPtr;

    // compute number of knots we can vary.  if this doesn't match the 
    // dimensions of variableKnots, we are screwed
    int nVarKnots=cls->ndata_sou-(cls->spline_degree+1);

    // set the knots 
    vector<double> myKnots(cls->ndata_sou+cls->spline_degree+1);
    for (int i=0;i<cls->spline_degree+1;++i) {
        myKnots[i]=cls->rmin;
        myKnots[cls->ndata_sou+cls->spline_degree-i]=cls->rmax;
    }
    for (int i=0;i<nVarKnots;++i) {
        while( (myKnots[cls->spline_degree+1+i]>cls->rmax) || (myKnots[cls->spline_degree+1+i]<cls->rmin) ) {
            myKnots[cls->spline_degree+1+i] =
                mod( abs(gsl_vector_get(variableKnots,i)), cls->rmax-cls->rmin ) + cls->rmin;
        }
    }
    sort(myKnots.begin(),myKnots.end());
    
    // Redimension souWork to accomodate the new knots
    cls->souWork.setDim( myKnots.size() - cls->souWork.spline_degree - 1 );

    // Now put in the new knots
    cls->souWork.knots = stl2tntVec(myKnots);

    // re-build the Kmtx
    parameterMap dummyMap;
    cls->set_kmtx(cls->corrwork,cls->souWork,dummyMap);
    
    // Now image the source.  
    return cls->imageit(cls->souWork); 

#else

    return 0.;

#endif

}

//------------------ set_knots --------------------------
//! master knot initializer routine
void CBasisSplineImager1d::set_knots(CBasisSpline1d& souout){

    // figure out if we have to use the default scheme instead of alternate
    if (knot_init_scheme=="sampling_thm") {
        if (spline_degree == 0){
            MESSAGE << "Requesting Spline degree = 0, so using default knots\n" <<ENDM_WARN;
            knot_init_scheme="default";
            ndata_sou=std::max(std::min(ndata_corr,ndata_sou),spline_degree+1);
        } 
        else if (kernel_particle1!=kernel_particle2){
            MESSAGE << "Sampling theorem doen't work on unlike pairs" <<ENDM_WARN;
            knot_init_scheme="default";
        } 
    }

    // these two make knots directly
    if (knot_init_scheme=="default") 
        set_knots_default(souout);
    else if (knot_init_scheme=="user_defined_knots") 
        set_knots_user_defined(souout);
    // these two make collocation points which must be converted to knots
    else if (knot_init_scheme=="sampling_thm") {
        set_colloc_sampling_thm();
        set_knots_from_colloc_points(souout);
    }
    else if (knot_init_scheme=="user_defined_collocation_points") 
        set_knots_from_colloc_points(souout);
    // default to this one
    else 
        set_knots_default(souout);
}

//------------------- set_knots_default -------------------------
void CBasisSplineImager1d::set_knots_default( CBasisSpline1d& souout ){

    cout << "  Default source knot initialization\n";
    // make a temp spline w/ correct dimensions, copy temp bspline into working source
    ndata_sou=std::min(static_cast<int>(ndata_corr/1.5),souout.ndata); //override if ndata_out unreasonable
    if (ndata_sou != souout.ndata){
        MESSAGE << "Number of coefficients requested is too big, lowering to "<<ndata_sou<<ENDM_WARN;
        souout.setDim( ndata_sou );
    }
    rmin = souout.xmin;
    rmax = souout.xmax;
    souout.setDefaultKnots();
}

//------------------- set_colloc_sampling_thm -------------------------
void CBasisSplineImager1d::set_colloc_sampling_thm( void ){

    static int maxnumzeros=40;

    cout << "   Using Sampling theorem to guess collocation points\n";
    cout << "       Getting zeroes from kernel\n";
    
    // get true zeros of wavefunction for collocation points
    vector<double> zeros(get_zeros(qscale,maxnumzeros,l));

    // make temp list ot r collocation points
    vector<double> tmp_rcolloc(1,rmin);

    // prune list of zeros to get only collocation points w/ r < _rmax
    for (int i=0;i<static_cast<int>(zeros.size());i++) {
        if ((zeros[i]<rmax)&&(zeros[i]>rmin)) tmp_rcolloc.push_back(zeros[i]);
    }

    // add last element
    tmp_rcolloc.push_back(rmax);

    // make sure list sorted
    sort(tmp_rcolloc.begin(),tmp_rcolloc.end());

    // save the collocation points
    colloc_pts = tmp_rcolloc;
}

//------------------- set_knots_from_colloc_points -------------------------
void CBasisSplineImager1d::set_knots_from_colloc_points( CBasisSpline1d& souout ){

    cout << endl;
    cout << "    Converting collocation points into knots\n";
    cout << "    Pre-thin collocation points:\n";
    for (int i=0;i<static_cast<int>(colloc_pts.size());++i){
        cout <<"        "<< i << " " << colloc_pts[i]<<endl;
    }

    // thin tmp_rcolloc
    vector<double> rcolloc(1,colloc_pts[0]);
    int nave=1;
    for (int i=1;i<static_cast<int>(colloc_pts.size());++i) {
        if (abs(colloc_pts[i]-rcolloc.back())<knot_tolerance) {
            rcolloc.back() = (nave/(nave+1))*(rcolloc.back()+colloc_pts[i]/nave);           
            nave++;
        } else {
            rcolloc.push_back(colloc_pts[i]);
            nave=1;
        }
    }
    cout << "    Post-thin collocation points:\n";
    for (int i=0;i<static_cast<int>(rcolloc.size());i++){
        cout << "        "<<i << " " << rcolloc[i]<<endl;
    }
    
    // create a temp bspline and set the optimal knots using the collocation points above
    souout.setOptimalKnots(stl2tntVec(rcolloc));
    // Since bspline knot setup may override the number of coefficients 
    // (we may have asked for unreasonable number), we need to reset ndata_out
    ndata_sou = souout.ndata;
}

//------------------- set_knots_user_defined -------------------------
void CBasisSplineImager1d::set_knots_user_defined( CBasisSpline1d& souout ){

    // copy temp bspline into working source
    souout.setDim(user_knots.size()-souout.spline_degree-1);
    // set the knots 
    souout.knots = stl2tntVec(user_knots);
    // Since bspline knot setup may override the number of coefficients 
    // (the user may have botched to knots/coeffs count), 
    // we need to reset ndata_out
    ndata_sou = souout.ndata;
}

//------------------- set_no_data -------------------------
bool CBasisSplineImager1d::set_no_data( CSourceFtnBase& souout ){

    CBasisFuncImager1d::set_no_data( souout );
    CBasisSpline1d* souPtr = dynamic_cast<CBasisSpline1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisSpline1d"<<ENDM_FATAL;
    souPtr->setDefaultKnots(); // will give 0 everywhere outside this range anyway
    return true;
}

//------------------- initialize_source -------------------------
bool CBasisSplineImager1d::initialize_source( const CCorrFtn1dHisto& corrin, CSourceFtnBase& souout, const parameterMap& m ){

    bool result = CBasisFuncImager1d::initialize_source(corrin,souout,m);
    CBasisSpline1d* souPtr = dynamic_cast<CBasisSpline1d*>(&souout);
    if (souPtr==NULL) throw MESSAGE<<"Source argument must derive from CBasisSpline1d"<<ENDM_FATAL;
    souPtr->Read(m);
    spline_degree=souPtr->spline_degree;
    set_knots(*souPtr);
    return result;
}

//----------------- get_zeros ---------------------------
//!  Get list of first nZeros zeros for l<=lmax=rmax*q/hbarc
vector< double > CBasisSplineImager1d::get_zeros(double q, int nZeros, int l, double eps){

    vector<double> zeros;
    
    // if I did this right, this should be smaller than the spacing between any zeros
    double xstep=0.5/(2.0*q/HBARC);  
    
    double xmin = eps;
    double testlo;
    double testhi;
    
    zero_struct azero;
    
    while (static_cast<int>(zeros.size())<nZeros){
        testlo  = kernelPtr->GetValue(l,q,xmin); 
        testhi  = kernelPtr->GetValue(l,q,xmin+xstep);

        // Make sure lower point isn't a zero
        if ( abs(testlo) <= eps ) { 
            zeros.push_back(xmin);
            xmin+=xstep/2.0;
            
        } else {

            // Make sure upper point isn't a zero
            if ( abs(testhi) <= eps ) { 
                testlo  = kernelPtr->GetValue(l,q,xmin);
                testhi  = kernelPtr->GetValue(l,q,xmin+xstep/2.0);
            }

            // OK, now look for zeros on the current interval
            if ( testlo*testhi <= 0.) {
                azero=find_zero(xmin,xmin+xstep,q,l,eps);
                if (azero.iszero) zeros.push_back(azero.result);
            }
        }
        
        // ready for next step just past this zero
        xmin+=xstep;
    }
    return zeros;
}

//------------------- find_zero ------------------------
//!  Gets first zero between xmin & xmax, 
//!  starting at xmin.  l and q are fixed. 
zero_struct CBasisSplineImager1d::find_zero(double xmin, double xmax, double q, int l, double eps){
    int nsteps=10;
    double stepsize=(xmax-xmin)/nsteps,xleft,xright;
    zero_struct ans;
    for (int i=0;i<nsteps;++i){
        xleft=xmin+i*stepsize;
        xright=xleft+stepsize;
        if (kernelPtr->GetValue(l,q,xleft)*kernelPtr->GetValue(l,q,xright)<=0.){ // is sign change
            if (abs(xleft-xright)<eps) {
                ans.iszero = true;
                ans.result = (xleft+xright)/2.;
                return ans;
            }
            else return find_zero(xleft,xright,q,l,eps);
        }
    }
    
    // never get here, unless is no zero
    ans.result=rmax; 
    ans.iszero=false; 
    return ans; 
}

