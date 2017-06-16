/* d07hre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int d07hre_(integer *ndim, integer *wtleng, doublereal *w, 
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
    static doublereal lam0, lam1, lam2, lamp, ratio, twondm;

/* ***BEGIN PROLOGUE D07HRE */
/* ***KEYWORDS basic integration rule, degree 7 */
/* ***PURPOSE  To initialize a degree 7 basic rule, and null rules. */
/* ***AUTHOR   Alan Genz, Computer Science Department, Washington */
/*            State University, Pullman, WA 99163-1210 USA */
/* ***LAST MODIFICATION 88-05-31 */
/* ***DESCRIPTION  D07HRE initializes a degree 7 integration rule, */
/*            two degree 5 null rules, one degree 3 null rule and one */
/*            degree 1 null rule for the hypercube [-1,1]**NDIM. */

/*   ON ENTRY */

/*   NDIM   Integer. */
/*          Number of variables. */
/*   WTLENG Integer. */
/*          The number of weights in each of the rules. */
/*          WTLENG MUST be set equal to 6. */

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
/* ***END PROLOGUE D07HRE */

/*   Global variables */


/*   Local Variables */


/* ***FIRST EXECUTABLE STATEMENT D07HRE */


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
    rulpts[*wtleng - 1] = (doublereal) ((*ndim << 1) * (*ndim - 1));
    rulpts[1] = 1.;

/*     Compute squared generator parameters */

    lam0 = (float).4707;
    lamp = (float).5625;
    lam1 = 4 / (15 - 5 / lam0);
    ratio = (1 - lam1 / lam0) / 27;
    lam2 = (5 - lam1 * 7 - ratio * 35) / (7 - lam1 * 35 / 3 - ratio * 35 / 
	    lam0);

/*     Compute degree 7 rule weights */

/* Computing 3rd power */
    d__1 = lam0 * 3;
    w[31] = 1 / (d__1 * (d__1 * d__1)) / twondm;
/* Computing 2nd power */
    d__1 = lam1;
    w[26] = (1 - lam0 * 5 / 3) / ((lam1 - lam0) * 60 * (d__1 * d__1));
    w[16] = (1 - lam2 * 5 / 3 - twondm * 5 * w[31] * lam0 * (lam0 - lam2)) / (
	    lam1 * 10 * (lam1 - lam2)) - (*ndim - 1 << 1) * w[26];
    w[11] = (1 - lam1 * 5 / 3 - twondm * 5 * w[31] * lam0 * (lam0 - lam1)) / (
	    lam2 * 10 * (lam2 - lam1));

/*     Compute weights for 2 degree 5, 1 degree 3 and 1 degree 1 rules */

/* Computing 3rd power */
    d__1 = lam0;
    w[32] = 1 / (d__1 * (d__1 * d__1) * 36) / twondm;
/* Computing 2nd power */
    d__1 = lam0;
/* Computing 2nd power */
    d__2 = lam1;
    w[27] = (1 - twondm * 9 * w[32] * (d__1 * d__1)) / (d__2 * d__2 * 36);
    w[17] = (1 - lam2 * 5 / 3 - twondm * 5 * w[32] * lam0 * (lam0 - lam2)) / (
	    lam1 * 10 * (lam1 - lam2)) - (*ndim - 1 << 1) * w[27];
    w[12] = (1 - lam1 * 5 / 3 - twondm * 5 * w[32] * lam0 * (lam0 - lam1)) / (
	    lam2 * 10 * (lam2 - lam1));
/* Computing 3rd power */
    d__1 = lam0;
    w[33] = 5 / (d__1 * (d__1 * d__1) * 108) / twondm;
/* Computing 2nd power */
    d__1 = lam0;
/* Computing 2nd power */
    d__2 = lam1;
    w[28] = (1 - twondm * 9 * w[33] * (d__1 * d__1)) / (d__2 * d__2 * 36);
    w[18] = (1 - lamp * 5 / 3 - twondm * 5 * w[33] * lam0 * (lam0 - lamp)) / (
	    lam1 * 10 * (lam1 - lamp)) - (*ndim - 1 << 1) * w[28];
    w[23] = (1 - lam1 * 5 / 3 - twondm * 5 * w[33] * lam0 * (lam0 - lam1)) / (
	    lamp * 10 * (lamp - lam1));
/* Computing 3rd power */
    d__1 = lam0;
    w[34] = 1 / (d__1 * (d__1 * d__1) * 54) / twondm;
/* Computing 2nd power */
    d__1 = lam0;
/* Computing 2nd power */
    d__2 = lam1;
    w[29] = (1 - twondm * 18 * w[34] * (d__1 * d__1)) / (d__2 * d__2 * 72);
    w[19] = (1 - lam2 * 10 / 3 - twondm * 10 * w[34] * lam0 * (lam0 - lam2)) /
	     (lam1 * 20 * (lam1 - lam2)) - (*ndim - 1 << 1) * w[29];
    w[14] = (1 - lam1 * 10 / 3 - twondm * 10 * w[34] * lam0 * (lam0 - lam1)) /
	     (lam2 * 20 * (lam2 - lam1));

/*     Set generator values */

    lam0 = sqrt(lam0);
    lam1 = sqrt(lam1);
    lam2 = sqrt(lam2);
    lamp = sqrt(lamp);
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__ + *wtleng * g_dim1] = lam0;
/* L40: */
    }
    g[(*wtleng - 1) * g_dim1 + 1] = lam1;
    g[(*wtleng - 1) * g_dim1 + 2] = lam1;
    g[(*wtleng - 4) * g_dim1 + 1] = lam2;
    g[(*wtleng - 3) * g_dim1 + 1] = lam1;
    g[(*wtleng - 2) * g_dim1 + 1] = lamp;

/*     Compute final weight values. */
/*     The null rule weights are computed from differences between */
/*     the degree 7 rule weights and lower degree rule weights. */

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
    errcof[3] = 1.;
    errcof[4] = 5.;
    errcof[5] = (float).5;
    errcof[6] = (float).25;

/* ***END D07HRE */

    return 0;
} /* d07hre_ */

#ifdef __cplusplus
	}
#endif
