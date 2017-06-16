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
#ifndef CORAL_H_INCLUDE_
#define CORAL_H_INCLUDE_
// This includes all headers

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cstring>

#include "utilities.h"
#include "kernel.h"
#include "source2cf.h"
#include "sourcecalc.h"
#include "wavefunction.h"
#include "cfcalc.h"
#include "sfit.h"

#include "kernel_chooser.h"
#include "convolution.h"
#include "expander.h"
#include "imager.h"
#include "general_imager3d.h"
#include "uncoupled_imager3d.h"
#include "general_imager1d.h"
#include "basisfunc_imager1d.h"
#include "bspline_imager1d.h"
#include "objects3d.h"
#include "histogram3d.h"
#include "harmonic_expansion.h"
#include "dataset.h"
#include "objects1d.h"
#include "func_expansion1d.h"
#include "generic_spline1d.h"
#include "histogram1d.h"
#include "basis_spline1d.h"
#include "orthogfunc_expansion1d.h"
#include "chebyshevpoly_expansion1d.h"
#include "hermitefunction_expansion1d.h"
#include "laguerrepoly_expansion1d.h"
#include "legendrepoly_expansion1d.h"
#include "corrbase.h"
#include "corr3d_histo.h"
#include "corr3d_ylm.h"
#include "corr1d_histo.h"
#include "soubase.h"
#include "sou3d_ylm.h"
#include "sou3d_histo.h"
#include "sou1d_gauss.h"
#include "sou1d_2gauss.h"
#include "sou1d_bsplines.h"
#include "sou1d_chebyshev.h"
#include "sou1d_histo.h"
#include "sou1d_hermite.h"
#include "sou1d_laguerre.h"
#include "sou1d_legendre.h"

#endif
