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
//
// We switched to using the cubature routine so this test is obsolete
//
// DAB 11/5/2009 
//




/* dtest1.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__5 = 5;
static integer c__6 = 6;
static integer c__5000 = 5000;
static integer c__0 = 0;
static integer c__1 = 1;


/*   DTEST1 is a simple test driver for DCUHRE. */

/*   Output produced on a SUN 3/50. */

/*       DCUHRE TEST RESULTS */

/*    FTEST CALLS = 3549, IFAIL =  0 */
/*   N   ESTIMATED ERROR   INTEGRAL */
/*   1     0.00000010     0.13850818 */
/*   2     0.00000013     0.06369469 */
/*   3     0.00000874     0.05861748 */
/*   4     0.00000021     0.05407034 */
/*   5     0.00000019     0.05005614 */
/*   6     0.00000009     0.04654608 */

/* Main program */ int MAIN__()
{
    /* Format strings */
    static char fmt_9999[] = "(8x,\002DCUHRE TEST RESULTS\002,//\002     FTE\
ST CALLS = \002,i4,\002, IFAIL = \002,i2,/\002    N   ESTIMATED ERROR   INTE\
GRAL\002)";
    static char fmt_9998[] = "(3x,i2,2f18.8)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static doublereal a[5], b[5];
    static integer n, key, ifail, neval;
    extern /* Subroutine */ int ftest_(...);
    extern /* Subroutine */ int dcuhre_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, U_fp, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *);
    static doublereal absreq, absest[6];
    static integer mincls, maxcls;
    static doublereal finest[6], relreq, wrkstr[5000];

    /* Fortran I/O blocks */
    static cilist io___14 = { 0, 6, 0, fmt_9999, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_9998, 0 };


    for (n = 1; n <= 5; ++n) {
	a[n - 1] = 0.;
	b[n - 1] = 1.;
/* L10: */
    }
    mincls = 0;
    maxcls = 10000;
    key = 0;
    absreq = 0.;
    relreq = (float).001;
    dcuhre_(&c__5, &c__6, a, b, &mincls, &maxcls, (U_fp)ftest_, &absreq, &
	    relreq, &key, &c__5000, &c__0, finest, absest, &neval, &ifail, 
	    wrkstr);
    s_wsfe(&io___14);
    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ifail, (ftnlen)sizeof(integer));
    e_wsfe();
    for (n = 1; n <= 6; ++n) {
	s_wsfe(&io___15);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&absest[n - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&finest[n - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L20: */
    }
    return 0;
} /* MAIN__ */

/* Subroutine */ int ftest_(integer *ndim, doublereal *z__, integer *nfun, 
	doublereal *f)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    double exp(doublereal);

    /* Local variables */
    static integer n;
    static doublereal sum;

    /* Parameter adjustments */
    --z__;
    --f;

    /* Function Body */
    sum = 0.;
    i__1 = *ndim;
    for (n = 1; n <= i__1; ++n) {
/* Computing 2nd power */
	d__1 = z__[n];
	sum += n * (d__1 * d__1);
/* L10: */
    }
    f[1] = exp(-sum / 2);
    i__1 = *ndim;
    for (n = 1; n <= i__1; ++n) {
	f[n + 1] = z__[n] * f[1];
/* L20: */
    }
    return 0;
} /* ftest_ */

/* Main program alias */ int dtest1_ () { MAIN__ (); return 0; }
#ifdef __cplusplus
	}
#endif
