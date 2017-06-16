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
#include "coral.h"

int main( void ){
    parameterMap m_in;
    parameterMap m_sou_in;
    parameterMap m_conv_in;
    parameterMap m_corr_in;
    parameterMap m_kern_in;
    // set up input source
    parameter::set(m_sou_in,"particle1",string("pi+"));
    parameter::set(m_sou_in,"particle2",string("pi+"));
    parameter::set(m_sou_in,"width",5.0);
    parameter::set(m_sou_in,"height",1.0);
    parameter::set(m_sou_in,"l",0);
    parameter::set(m_sou_in,"m",0);
    parameter::set(m_sou_in,"realpart",true);
    // set up for the convoluter
    parameter::set(m_conv_in,"rmin",0.0);
    parameter::set(m_conv_in,"rmax",50.0);
    // set up for the kernel
    parameter::set(m_kern_in,"hbt_only",false);
    parameter::set(m_kern_in,"kernel_cache",string("kcache/pi+pi+/"));
    parameter::set(m_kern_in,"use_cache",true);
    parameter::set(m_kern_in,"read_cache",true);
    parameter::set(m_kern_in,"param_filename",string("wfparameters.dat"));
    // set up for the correlation we're creating
    parameter::set(m_corr_in,"fixed_width_bins",true);
    parameter::set(m_corr_in,"dx",2.55);
    parameter::set(m_corr_in,"xoffset",1.2750);
    parameter::set(m_corr_in,"ndata",40);
    // load 'em into the main map
    parameter::set(m_in,"convoluter_settings",m_conv_in);
    parameter::set(m_in,"correlation_settings",m_corr_in);
    parameter::set(m_in,"kernel_settings",m_kern_in);

    CGaussianSource source;
    source.Read(m_sou_in);
    parameterMap m_test;
    source.Write(m_test);
    cout << "_____Gaussian Source Parameters_____";
    cout << m_test <<endl;
    cout << "source.getValue(2.0)="<<source.getValue(2.0)<<endl <<endl;
    
    CCorrFtn1dHisto correlation(convolute(source,m_in));
    for (int i=0;i<correlation.ndata;++i){
        correlation.uncert[i]=correlation.data[i]*0.1;
    }
    
    parameterMap m_out;
    correlation.Write(m_out);
    ofstream corrout("pion_correlation.dat");
    corrout << m_out;

    cout << "\nSample points from correlation:\nq   C(q)"<<endl;
    for (int i=0;i<10;++i){
        double q=5.0*(double(i)+0.5);
        cout << q << " " << correlation.getValue(q) << endl;
    }
    
    return true;
}
