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





/* dtest2.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"
#include <stdio.h>

/* Table of constant values */

static integer c__2 = 2;
static integer c__195 = 195;
static integer c__16 = 16;
static integer c__1 = 1;
static integer c__5000 = 5000;
static integer c__0 = 0;
static integer c__3 = 3;
static integer c__381 = 381;
static integer c__4 = 4;
static integer c__459 = 459;
static integer c__5 = 5;
static integer c__309 = 309;
static integer c__6 = 6;
static integer c__3000 = 3000;
static integer c__6000 = 6000;
static integer c__12000 = 12000;
static integer c__500 = 500;
static integer c__1500 = 1500;
static integer c__7 = 7;
static integer c__10000 = 10000;
static integer c__10 = 10;
static integer c__20000 = 20000;

/* Main program */ int MAIN__()
{
    static doublereal a[10], b[10];
    static integer n;
    extern /* Subroutine */ int atest_(integer *, doublereal *, doublereal *, 
	    integer *, integer *, U_fp, doublereal *, integer *, integer *, 
	    integer *, doublereal *);
    static doublereal abserr;
    extern /* Subroutine */ int ftesto_(...), ftestp_(...), ftestx_(...);
    static doublereal wrkstr[5000];


/*   DTEST2 tests some of the features of DCUHRE. */
/*   DTEST2 checks that DCUHRE integrates to machine */
/*   precision some of the monomials that DCUHRE is */
/*   supposed to integrate to machine precision. */
/*   DTEST2 checks that the restart feature of DCUHRE works. */
/*   DTEST2 runs small tests in dimensions 2, 3, 5, 7 and 10. */

/*   Output produced on a SUN 3/50. */



/*    DCUHRE TEST WITH NDIM =   2, KEY =  1 */
/*    SUBROUTINE CALLS =    195, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.00e+00   1.0000000000000000 */
/*      2      0.22e-15   1.0000000000000002 */
/*      3      0.22e-15   1.0000000000000002 */
/*      4      0.22e-15   1.0000000000000002 */
/*      5      0.22e-15   1.0000000000000002 */
/*      6      0.00e+00   1.0000000000000000 */
/*      7      0.22e-15   1.0000000000000002 */
/*      8      0.00e+00   1.0000000000000000 */
/*      9      0.00e+00   1.0000000000000000 */
/*     10      0.22e-15   1.0000000000000002 */
/*     11      0.00e+00   1.0000000000000000 */
/*     12      0.22e-15   1.0000000000000002 */
/*     13      0.22e-15   1.0000000000000002 */
/*     14      0.22e-15   1.0000000000000002 */
/*     15      0.22e-15   1.0000000000000002 */
/*     16      0.22e-15   1.0000000000000002 */

/*    DCUHRE TEST WITH NDIM =   3, KEY =  2 */
/*    SUBROUTINE CALLS =    381, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.00e+00   1.0000000000000000 */
/*      2      0.00e+00   1.0000000000000000 */
/*      3      0.00e+00   1.0000000000000000 */
/*      4      0.00e+00   1.0000000000000000 */
/*      5      0.11e-15   0.9999999999999999 */
/*      6      0.00e+00   1.0000000000000000 */
/*      7      0.00e+00   1.0000000000000000 */
/*      8      0.11e-15   0.9999999999999999 */
/*      9      0.00e+00   1.0000000000000000 */
/*     10      0.00e+00   1.0000000000000000 */
/*     11      0.00e+00   1.0000000000000000 */
/*     12      0.00e+00   1.0000000000000000 */
/*     13      0.11e-15   0.9999999999999999 */
/*     14      0.00e+00   1.0000000000000000 */
/*     15      0.00e+00   1.0000000000000000 */
/*     16      0.87e-09   1.0000000008661680 */

/*    DCUHRE TEST WITH NDIM =   4, KEY =  3 */
/*    SUBROUTINE CALLS =    459, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.44e-15   1.0000000000000004 */
/*      2      0.00e+00   1.0000000000000000 */
/*      3      0.44e-15   1.0000000000000004 */
/*      4      0.00e+00   1.0000000000000000 */
/*      5      0.00e+00   1.0000000000000000 */
/*      6      0.44e-15   1.0000000000000004 */
/*      7      0.22e-15   1.0000000000000002 */
/*      8      0.22e-15   1.0000000000000002 */
/*      9      0.11e-15   0.9999999999999999 */
/*     10      0.22e-15   1.0000000000000002 */
/*     11      0.44e-15   1.0000000000000004 */
/*     12      0.44e-15   1.0000000000000004 */
/*     13      0.00e+00   1.0000000000000000 */
/*     14      0.11e-15   0.9999999999999999 */
/*     15      0.24e-07   0.9999999758753315 */
/*     16      0.33e-06   0.9999996739089504 */

/*    DCUHRE TEST WITH NDIM =   5, KEY =  4 */
/*    SUBROUTINE CALLS =    309, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.00e+00   1.0000000000000000 */
/*      2      0.00e+00   1.0000000000000000 */
/*      3      0.00e+00   1.0000000000000000 */
/*      4      0.00e+00   1.0000000000000000 */
/*      5      0.22e-15   0.9999999999999998 */
/*      6      0.00e+00   1.0000000000000000 */
/*      7      0.11e-15   0.9999999999999999 */
/*      8      0.00e+00   1.0000000000000000 */
/*      9      0.22e-15   1.0000000000000002 */
/*     10      0.22e-15   1.0000000000000002 */
/*     11      0.22e-15   1.0000000000000002 */
/*     12      0.22e-15   1.0000000000000002 */
/*     13      0.22e-15   1.0000000000000002 */
/*     14      0.15e-05   1.0000015013499370 */
/*     15      0.15e-04   1.0000147498974314 */
/*     16      0.76e-04   1.0000756306700240 */

/*    DCUHRE TEST WITH NDIM =   6, KEY =  4 */
/*    SUBROUTINE CALLS =   2737, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.23e-05   1.0000023460841487 */

/*    DCUHRE TEST WITH NDIM =   6, KEY =  4 */
/*    SUBROUTINE CALLS =   5957, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.11e-05   1.0000010801329045 */

/*    DCUHRE TEST WITH NDIM =   6, KEY =  4 */
/*    SUBROUTINE CALLS =  11753, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.62e-06   1.0000006230062635 */

/*    DCUHRE TEST WITH NDIM =   2, KEY =  1 */
/*    SUBROUTINE CALLS =    455, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.19e-09   0.9999999998085931 */

/*    DCUHRE TEST WITH NDIM =   3, KEY =  2 */
/*    SUBROUTINE CALLS =   1397, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.57e-07   1.0000000565042295 */

/*    DCUHRE TEST WITH NDIM =   5, KEY =  3 */
/*    SUBROUTINE CALLS =   4641, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.13e-05   0.9999987098089600 */

/*    DCUHRE TEST WITH NDIM =   7, KEY =  4 */
/*    SUBROUTINE CALLS =   9945, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.35e-06   1.0000003517922749 */

/*    DCUHRE TEST WITH NDIM =  10, KEY =  4 */
/*    SUBROUTINE CALLS =  18975, IFAIL =  1 */
/*      N   ABSOLUTE ERROR    INTEGRAL */
/*      1      0.39e-04   0.9999614789438658 */

    for (n = 1; n <= 10; ++n) {
	a[n - 1] = 0.;
	b[n - 1] = 1.;
/* L10: */
    }
    abserr = (float)1e-10;

/*    TEST FOR INTEGRATING POLYNOMIALS */
/*     Selected monomials, degrees 0-13 */

/*         Degree 13 rule */
    atest_(&c__2, a, b, &c__195, &c__16, (U_fp)ftestp_, &abserr, &c__1, &
	    c__5000, &c__0, wrkstr);
/*         Degree 11 rule */
    atest_(&c__3, a, b, &c__381, &c__16, (U_fp)ftestp_, &abserr, &c__2, &
	    c__5000, &c__0, wrkstr);
/*         Degree  9 rule */
    atest_(&c__4, a, b, &c__459, &c__16, (U_fp)ftestp_, &abserr, &c__3, &
	    c__5000, &c__0, wrkstr);
/*         Degree  7 rule */
    atest_(&c__5, a, b, &c__309, &c__16, (U_fp)ftestp_, &abserr, &c__4, &
	    c__5000, &c__0, wrkstr);

/*    TEST RESTART */

    atest_(&c__6, a, b, &c__3000, &c__1, (U_fp)ftesto_, &abserr, &c__4, &
	    c__5000, &c__0, wrkstr);
    atest_(&c__6, a, b, &c__6000, &c__1, (U_fp)ftesto_, &abserr, &c__4, &
	    c__5000, &c__1, wrkstr);
    atest_(&c__6, a, b, &c__12000, &c__1, (U_fp)ftesto_, &abserr, &c__4, &
	    c__5000, &c__1, wrkstr);

/*    TEST WITH NDIM = 2, 3, 5, 7, 10 */

    atest_(&c__2, a, b, &c__500, &c__1, (U_fp)ftestx_, &abserr, &c__1, &
	    c__5000, &c__0, wrkstr);
    atest_(&c__3, a, b, &c__1500, &c__1, (U_fp)ftestx_, &abserr, &c__2, &
	    c__5000, &c__0, wrkstr);
    atest_(&c__5, a, b, &c__5000, &c__1, (U_fp)ftestx_, &abserr, &c__3, &
	    c__5000, &c__0, wrkstr);
    atest_(&c__7, a, b, &c__10000, &c__1, (U_fp)ftestx_, &abserr, &c__4, &
	    c__5000, &c__0, wrkstr);
    atest_(&c__10, a, b, &c__20000, &c__1, (U_fp)ftestx_, &abserr, &c__4, &
	    c__5000, &c__0, wrkstr);

    return 0;
} /* MAIN__ */

/* Subroutine */ int atest_(integer *ndim, doublereal *a, doublereal *b, 
	integer *maxcls, integer *nfun, U_fp tstsub, doublereal *abserr, 
	integer *key, integer *lenwrk, integer *irest, doublereal *wrkstr)
{
    /* Format strings */
    static char fmt_99999[] = "(/5x,\002DCUHRE TEST WITH NDIM = \002,i3,\002\
, KEY = \002,i2)";
    static char fmt_99998[] = "(5x,\002SUBROUTINE CALLS = \002,i12,\002, IFAI\
L = \002,i2)";
    static char fmt_99997[] = "(7x,\002N   ABSOLUTE ERROR    INTEGRAL\002)";
    static char fmt_99996[] = "(6x,i2,e16.2,2x,f28.16)";

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();

    /* Local variables */
    static integer n;
    static doublereal rel;
    static integer ifail, neval;
    extern /* Subroutine */ int dcuhre_(integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, U_fp, doublereal *, 
	    doublereal *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, doublereal *);
    static doublereal absest[20], finest[20];

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, fmt_99999, 0 };
    static cilist io___12 = { 0, 6, 0, fmt_99998, 0 };
    static cilist io___13 = { 0, 6, 0, fmt_99997, 0 };
    static cilist io___15 = { 0, 6, 0, fmt_99996, 0 };


    /* Parameter adjustments */
    --b;
    --a;
    --wrkstr;

    /* Function Body */
    rel = 0.;
    dcuhre_(ndim, nfun, &a[1], &b[1], &c__0, maxcls, (U_fp)tstsub, abserr, &
	    rel, key, lenwrk, irest, finest, absest, &neval, &ifail, &wrkstr[
	    1]);

printf( "\n" );
printf( "ndim %i", *ndim );
printf( ", nfun %i", *nfun );
printf( ", neval %i", neval );
printf( ", maxcls %i", *maxcls );
printf( ", ifail %i", ifail );
printf( ", key %i", *key );
printf( ", lenwrk %i", *lenwrk );
printf( ", irest %i", *irest );
printf( ", wrkstr %e", *wrkstr );
printf( ", a %e", *a );
printf( ", b %e", *b );
printf( ", abserr %e", *abserr );
printf( "\n" );

    s_wsfe(&io___11);
    do_fio(&c__1, (char *)&(*ndim), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*key), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___12);
    do_fio(&c__1, (char *)&neval, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&ifail, (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___13);
    e_wsfe();
    i__1 = *nfun;
    for (n = 1; n <= i__1; ++n) {
	s_wsfe(&io___15);
	do_fio(&c__1, (char *)&n, (ftnlen)sizeof(integer));
	d__2 = (d__1 = finest[n - 1] - 1, abs(d__1));
	do_fio(&c__1, (char *)&d__2, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&finest[n - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
/* L10: */
    }
    return 0;
} /* atest_ */

/* Subroutine */ int ftestp_(integer *ndim, doublereal *z__, integer *nfun, 
	doublereal *f)
{
    /* System generated locals */
    doublereal d__1, d__2;


/*       Selected monomials, degree 0-13 */

    /* Parameter adjustments */
    --z__;
    --f;

    /* Function Body */
    f[1] = 1.;
    f[2] = z__[1] * 2;
/* Computing 2nd power */
    d__1 = z__[1];
    f[3] = d__1 * d__1 * 3;
    f[4] = f[2] * 2 * z__[2];
/* Computing 3rd power */
    d__1 = z__[1];
    f[5] = d__1 * (d__1 * d__1) * 4;
    f[6] = f[3] * 2 * z__[2];
/* Computing 4th power */
    d__1 = z__[1], d__1 *= d__1;
    f[7] = d__1 * d__1 * 5;
    f[8] = f[5] * 2 * z__[2];
/* Computing 2nd power */
    d__1 = z__[2];
    f[9] = f[3] * 3 * (d__1 * d__1);
/* Computing 5th power */
    d__1 = z__[1], d__2 = d__1, d__1 *= d__1;
    f[10] = d__2 * (d__1 * d__1) * 6;
    f[11] = f[7] * 2 * z__[2];
/* Computing 2nd power */
    d__1 = z__[2];
    f[12] = f[5] * 3 * (d__1 * d__1);
/* Computing 7th power */
    d__1 = z__[1], d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
    f[13] = d__2 * (d__1 * d__1) * 8;
/* Computing 9th power */
    d__1 = z__[1], d__2 = d__1, d__1 *= d__1, d__1 *= d__1;
    f[14] = d__2 * (d__1 * d__1) * 10;
/* Computing 11th power */
    d__1 = z__[1], d__2 = d__1, d__1 *= d__1, d__2 *= d__1, d__1 *= d__1;
    f[15] = d__2 * (d__1 * d__1) * 12;
/* Computing 13th power */
    d__1 = z__[1], d__2 = d__1, d__1 *= d__1, d__1 *= d__1, d__2 *= d__1;
    f[16] = d__2 * (d__1 * d__1) * 14;
    return 0;
} /* ftestp_ */

/* Subroutine */ int ftesto_(integer *ndim, doublereal *z__, integer *nfun, 
	doublereal *f)
{
    /* System generated locals */
    doublereal d__1;


/*     Corner Peak */

    /* Parameter adjustments */
    --z__;
    --f;

    /* Function Body */
/* Computing 6th power */
    d__1 = z__[1] * (float).1 + 1 + z__[2] * (float).2 + z__[3] * (float).3 + 
	    z__[4] * (float).4 + z__[5] * (float).5 + z__[6] * (float).6, 
	    d__1 *= d__1;
    f[1] = 10 / (d__1 * (d__1 * d__1)) / (float).2057746;
    return 0;
} /* ftesto_ */

/* Subroutine */ int ftestx_(integer *ndim, doublereal *z__, integer *nfun, 
	doublereal *f)
{
    /* Builtin functions */
    double cos(doublereal);

    /* Local variables */
    static integer n;
    static doublereal sum;


/*     Sum of Cosines */

    /* Parameter adjustments */
    --z__;
    --f;

    /* Function Body */
    sum = 0.;
    for (n = 1; n <= 2; ++n) {
	sum -= cos(z__[n] * 10) / (float).054402111088937;
/* L10: */
    }
    f[1] = sum / 2;
    return 0;
} /* ftestx_ */

/* Main program alias */ int dtest2_ () { MAIN__ (); return 0; }
#ifdef __cplusplus
	}
#endif
