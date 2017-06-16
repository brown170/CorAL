/* d09hre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int d09hre_(integer *ndim, integer *wtleng, doublereal *w, 
	doublereal *g, doublereal *errcof, doublereal *rulpts)
{
    /* System generated locals */
    integer g_dim1, g_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    //integer pow_ii(integer *, integer *);
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal lam0, lam1, lam2, lam3, lamp, ratio, twondm;

/* ***BEGIN PROLOGUE D09HRE */
/* ***KEYWORDS basic integration rule, degree 9 */
/* ***PURPOSE  To initialize a degree 9 basic rule and null rules. */
/* ***AUTHOR   Alan Genz, Computer Science Department, Washington */
/*            State University, Pullman, WA 99163-1210 USA */
/* ***LAST MODIFICATION 88-05-20 */
/* ***DESCRIPTION  D09HRE initializes a degree 9 integration rule, */
/*            two degree 7 null rules, one degree 5 null rule and one */
/*            degree 3 null rule for the hypercube [-1,1]**NDIM. */

/*   ON ENTRY */

/*   NDIM   Integer. */
/*          Number of variables. */
/*   WTLENG Integer. */
/*          The number of weights in each of the rules. */

/*   ON RETURN */
/*   W      Real array of dimension (5,WTLENG). */
/*          The weights for the basic and null rules. */
/*          W(1,1),...,W(1,WTLENG) are weights for the basic rule. */
/*          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights. */
/*   G      Real array of dimension (NDIM, WTLENG). */
/*          The fully symmetric sum generators for the rules. */
/*          G(1, J), ..., G(NDIM, J) are the are the generators for the */
/*          points associated with the Jth weights. */
/*   ERRCOF Real array of dimension 6. */
/*          Heuristic error coefficients that are used in the */
/*          error estimation in BASRUL. */
/*   RULPTS Real array of dimension WTLENG. */
/*          A work array. */

/* ***REFERENCES A. Genz and A. Malik, */
/*             "An Imbedded Family of Fully Symmetric Numerical */
/*              Integration Rules", */
/*              SIAM J Numer. Anal. 20 (1983), pp. 580-588. */
/* ***ROUTINES CALLED-NONE */
/* ***END PROLOGUE D09HRE */

/*   Global variables */


/*   Local Variables */


/* ***FIRST EXECUTABLE STATEMENT D09HRE */


/*     Initialize generators, weights and RULPTS */

    /* Parameter adjustments */
    --rulpts;
    g_dim1 = *ndim;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    w -= 6;
    --errcof;

    /* Function Body */
    i__1 = *wtleng;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *ndim;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    g[i__ + j * g_dim1] = 0.;
/* L10: */
	}
	for (i__ = 1; i__ <= 5; ++i__) {
	    w[i__ + j * 5] = 0.;
/* L20: */
	}
	rulpts[j] = (doublereal) (*ndim << 1);
/* L30: */
    }
    //twondm = (doublereal) pow_ii(&c__2, ndim);
    twondm = (doublereal) pow_ii(c__2, *ndim);
    rulpts[*wtleng] = twondm;
    if (*ndim > 2) {
	rulpts[8] = (doublereal) ((*ndim << 2) * (*ndim - 1) * (*ndim - 2) / 
		3);
    }
    rulpts[7] = (doublereal) ((*ndim << 2) * (*ndim - 1));
    rulpts[6] = (doublereal) ((*ndim << 1) * (*ndim - 1));
    rulpts[1] = 1.;

/*     Compute squared generator parameters */

    lam0 = (float).4707;
    lam1 = 4 / (15 - 5 / lam0);
    ratio = (1 - lam1 / lam0) / 27;
    lam2 = (5 - lam1 * 7 - ratio * 35) / (7 - lam1 * 35 / 3 - ratio * 35 / 
	    lam0);
    ratio = ratio * (1 - lam2 / lam0) / 3;
    lam3 = (7 - (lam2 + lam1) * 9 + lam2 * 63 * lam1 / 5 - ratio * 63) / (9 - 
	    (lam2 + lam1) * 63 / 5 + lam2 * 21 * lam1 - ratio * 63 / lam0);
    lamp = (float).0625;

/*     Compute degree 9 rule weights */

/* Computing 4th power */
    d__1 = lam0 * 3, d__1 *= d__1;
    w[*wtleng * 5 + 1] = 1 / (d__1 * d__1) / twondm;
    if (*ndim > 2) {
/* Computing 3rd power */
	d__1 = lam1 * 6;
	w[41] = (1 - 1 / (lam0 * 3)) / (d__1 * (d__1 * d__1));
    }
    w[36] = (1 - (lam0 + lam1) * 7 / 5 + lam0 * 7 * lam1 / 3) / (lam1 * 84 * 
	    lam2 * (lam2 - lam0) * (lam2 - lam1));
    w[31] = (1 - (lam0 + lam2) * 7 / 5 + lam0 * 7 * lam2 / 3) / (lam1 * 84 * 
	    lam1 * (lam1 - lam0) * (lam1 - lam2)) - w[36] * lam2 / lam1 - (*
	    ndim - 2 << 1) * w[41];
    w[21] = (1 - ((lam0 + lam1 + lam2) / 7 - (lam0 * lam1 + lam0 * lam2 + 
	    lam1 * lam2) / 5) * 9 - lam0 * 3 * lam1 * lam2) / (lam3 * 18 * (
	    lam3 - lam0) * (lam3 - lam1) * (lam3 - lam2));
    w[16] = (1 - ((lam0 + lam1 + lam3) / 7 - (lam0 * lam1 + lam0 * lam3 + 
	    lam1 * lam3) / 5) * 9 - lam0 * 3 * lam1 * lam3) / (lam2 * 18 * (
	    lam2 - lam0) * (lam2 - lam1) * (lam2 - lam3)) - (*ndim - 1 << 1) *
	     w[36];
    w[11] = (1 - ((lam0 + lam2 + lam3) / 7 - (lam0 * lam2 + lam0 * lam3 + 
	    lam2 * lam3) / 5) * 9 - lam0 * 3 * lam2 * lam3) / (lam1 * 18 * (
	    lam1 - lam0) * (lam1 - lam2) * (lam1 - lam3)) - (*ndim - 1 << 1) *
	     (w[36] + w[31] + (*ndim - 2) * w[41]);

/*     Compute weights for 2 degree 7, 1 degree 5 and 1 degree 3 rules */

/* Computing 4th power */
    d__1 = lam0, d__1 *= d__1;
    w[*wtleng * 5 + 2] = 1 / (d__1 * d__1 * 108) / twondm;
    if (*ndim > 2) {
/* Computing 3rd power */
	d__1 = lam0;
/* Computing 3rd power */
	d__2 = lam1 * 6;
	w[42] = (1 - twondm * 27 * w[47] * (d__1 * (d__1 * d__1))) / (d__2 * (
		d__2 * d__2));
    }
/* Computing 2nd power */
    d__1 = lam0;
    w[37] = (1 - lam1 * 5 / 3 - twondm * 15 * w[*wtleng * 5 + 2] * (d__1 * 
	    d__1) * (lam0 - lam1)) / (lam1 * 60 * lam2 * (lam2 - lam1));
/* Computing 2nd power */
    d__1 = lam0;
    w[32] = (1 - (lam1 * 8 * lam2 * w[37] + twondm * w[*wtleng * 5 + 2] * (
	    d__1 * d__1)) * 9) / (lam1 * 36 * lam1) - w[42] * 2 * (*ndim - 2);
    w[22] = (1 - ((lam1 + lam2) / 5 - lam1 * lam2 / 3 + twondm * w[*wtleng * 
	    5 + 2] * lam0 * (lam0 - lam1) * (lam0 - lam2)) * 7) / (lam3 * 14 *
	     (lam3 - lam1) * (lam3 - lam2));
    w[17] = (1 - ((lam1 + lam3) / 5 - lam1 * lam3 / 3 + twondm * w[*wtleng * 
	    5 + 2] * lam0 * (lam0 - lam1) * (lam0 - lam3)) * 7) / (lam2 * 14 *
	     (lam2 - lam1) * (lam2 - lam3)) - (*ndim - 1 << 1) * w[37];
    w[12] = (1 - ((lam2 + lam3) / 5 - lam2 * lam3 / 3 + twondm * w[*wtleng * 
	    5 + 2] * lam0 * (lam0 - lam2) * (lam0 - lam3)) * 7) / (lam1 * 14 *
	     (lam1 - lam2) * (lam1 - lam3)) - (*ndim - 1 << 1) * (w[37] + w[
	    32] + (*ndim - 2) * w[42]);
/* Computing 4th power */
    d__1 = lam0, d__1 *= d__1;
    w[*wtleng * 5 + 3] = 5 / (d__1 * d__1 * 324) / twondm;
    if (*ndim > 2) {
/* Computing 3rd power */
	d__1 = lam0;
/* Computing 3rd power */
	d__2 = lam1 * 6;
	w[43] = (1 - twondm * 27 * w[48] * (d__1 * (d__1 * d__1))) / (d__2 * (
		d__2 * d__2));
    }
/* Computing 2nd power */
    d__1 = lam0;
    w[38] = (1 - lam1 * 5 / 3 - twondm * 15 * w[*wtleng * 5 + 3] * (d__1 * 
	    d__1) * (lam0 - lam1)) / (lam1 * 60 * lam2 * (lam2 - lam1));
/* Computing 2nd power */
    d__1 = lam0;
    w[33] = (1 - (lam1 * 8 * lam2 * w[38] + twondm * w[*wtleng * 5 + 3] * (
	    d__1 * d__1)) * 9) / (lam1 * 36 * lam1) - w[43] * 2 * (*ndim - 2);
    w[28] = (1 - ((lam1 + lam2) / 5 - lam1 * lam2 / 3 + twondm * w[*wtleng * 
	    5 + 3] * lam0 * (lam0 - lam1) * (lam0 - lam2)) * 7) / (lamp * 14 *
	     (lamp - lam1) * (lamp - lam2));
    w[18] = (1 - ((lam1 + lamp) / 5 - lam1 * lamp / 3 + twondm * w[*wtleng * 
	    5 + 3] * lam0 * (lam0 - lam1) * (lam0 - lamp)) * 7) / (lam2 * 14 *
	     (lam2 - lam1) * (lam2 - lamp)) - (*ndim - 1 << 1) * w[38];
    w[13] = (1 - ((lam2 + lamp) / 5 - lam2 * lamp / 3 + twondm * w[*wtleng * 
	    5 + 3] * lam0 * (lam0 - lam2) * (lam0 - lamp)) * 7) / (lam1 * 14 *
	     (lam1 - lam2) * (lam1 - lamp)) - (*ndim - 1 << 1) * (w[38] + w[
	    33] + (*ndim - 2) * w[43]);
/* Computing 4th power */
    d__1 = lam0, d__1 *= d__1;
    w[*wtleng * 5 + 4] = 2 / (d__1 * d__1 * 81) / twondm;
    if (*ndim > 2) {
/* Computing 3rd power */
	d__1 = lam0;
/* Computing 3rd power */
	d__2 = lam1 * 6;
	w[44] = (2 - twondm * 27 * w[49] * (d__1 * (d__1 * d__1))) / (d__2 * (
		d__2 * d__2));
    }
    w[39] = (2 - lam1 * 15 / 9 - twondm * 15 * w[*wtleng * 5 + 4] * lam0 * (
	    lam0 - lam1)) / (lam1 * 60 * lam2 * (lam2 - lam1));
/* Computing 2nd power */
    d__1 = lam0;
    w[34] = (1 - (lam1 * 8 * lam2 * w[39] + twondm * w[*wtleng * 5 + 4] * (
	    d__1 * d__1)) * 9) / (lam1 * 36 * lam1) - w[44] * 2 * (*ndim - 2);
    w[24] = (2 - ((lam1 + lam2) / 5 - lam1 * lam2 / 3 + twondm * w[*wtleng * 
	    5 + 4] * lam0 * (lam0 - lam1) * (lam0 - lam2)) * 7) / (lam3 * 14 *
	     (lam3 - lam1) * (lam3 - lam2));
    w[19] = (2 - ((lam1 + lam3) / 5 - lam1 * lam3 / 3 + twondm * w[*wtleng * 
	    5 + 4] * lam0 * (lam0 - lam1) * (lam0 - lam3)) * 7) / (lam2 * 14 *
	     (lam2 - lam1) * (lam2 - lam3)) - (*ndim - 1 << 1) * w[39];
    w[14] = (2 - ((lam2 + lam3) / 5 - lam2 * lam3 / 3 + twondm * w[*wtleng * 
	    5 + 4] * lam0 * (lam0 - lam2) * (lam0 - lam3)) * 7) / (lam1 * 14 *
	     (lam1 - lam2) * (lam1 - lam3)) - (*ndim - 1 << 1) * (w[39] + w[
	    34] + (*ndim - 2) * w[44]);
    w[15] = 1 / (lam1 * 6);

/*     Set generator values */

    lam0 = sqrt(lam0);
    lam1 = sqrt(lam1);
    lam2 = sqrt(lam2);
    lam3 = sqrt(lam3);
    lamp = sqrt(lamp);
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__ + *wtleng * g_dim1] = lam0;
/* L40: */
    }
    if (*ndim > 2) {
	g[(g_dim1 << 3) + 1] = lam1;
	g[(g_dim1 << 3) + 2] = lam1;
	g[(g_dim1 << 3) + 3] = lam1;
    }
    g[g_dim1 * 7 + 1] = lam1;
    g[g_dim1 * 7 + 2] = lam2;
    g[g_dim1 * 6 + 1] = lam1;
    g[g_dim1 * 6 + 2] = lam1;
    g[g_dim1 * 5 + 1] = lamp;
    g[(g_dim1 << 2) + 1] = lam3;
    g[g_dim1 * 3 + 1] = lam2;
    g[(g_dim1 << 1) + 1] = lam1;

/*     Compute final weight values. */
/*     The null rule weights are computed from differences between */
/*     the degree 9 rule weights and lower degree rule weights. */

    w[6] = twondm;
    for (j = 2; j <= 5; ++j) {
	i__1 = *wtleng;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    w[j + i__ * 5] -= w[i__ * 5 + 1];
	    w[j + 5] -= rulpts[i__] * w[j + i__ * 5];
/* L50: */
	}
/* L70: */
    }
    i__1 = *wtleng;
    for (i__ = 2; i__ <= i__1; ++i__) {
	w[i__ * 5 + 1] = twondm * w[i__ * 5 + 1];
	w[6] -= rulpts[i__] * w[i__ * 5 + 1];
/* L80: */
    }

/*     Set error coefficients */

    errcof[1] = 5.;
    errcof[2] = 5.;
    errcof[3] = (float)1.;
    errcof[4] = 5.;
    errcof[5] = (float).5;
    errcof[6] = (float).25;

/* ***END D09HRE */

    return 0;
} /* d09hre_ */

#ifdef __cplusplus
	}
#endif
