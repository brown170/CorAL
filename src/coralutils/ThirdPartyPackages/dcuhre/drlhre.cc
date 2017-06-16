/* drlhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int drlhre_(integer *ndim, doublereal *center, doublereal *
	hwidth, integer *wtleng, doublereal *g, doublereal *w, doublereal *
	errcof, integer *numfun, S_fp funsub, doublereal *scales, doublereal *
	norms, doublereal *x, doublereal *null, doublereal *basval, 
	doublereal *rgnerr, doublereal *direct)
{
    /* System generated locals */
    integer null_dim1, null_offset, g_dim1, g_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Local variables */
    static integer i__, j, k;
    static doublereal ratio, search, difmax;
    extern /* Subroutine */ int dfshre_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, S_fp, doublereal *, 
	    doublereal *);
    static doublereal frthdf, difsum;
    static integer divaxn;
    static doublereal rgnvol;

/* ***BEGIN PROLOGUE DRLHRE */
/* ***KEYWORDS basic numerical integration rule */
/* ***PURPOSE  To compute basic integration rule values. */
/* ***AUTHOR   Alan Genz, Computer Science Department, Washington */
/*            State University, Pullman, WA 99163-1210 USA */
/* ***LAST MODIFICATION 90-02-06 */
/* ***DESCRIPTION DRLHRE computes basic integration rule values for a */
/*            vector of integrands over a hyper-rectangular region. */
/*            These are estimates for the integrals. DRLHRE also computes */
/*            estimates for the errors and determines the coordinate axis */
/*            where the fourth difference for the integrands is largest. */

/*   ON ENTRY */

/*   NDIM   Integer. */
/*          Number of variables. */
/*   CENTER Real array of dimension NDIM. */
/*          The coordinates for the center of the region. */
/*   HWIDTH Real Array of dimension NDIM. */
/*          HWIDTH(I) is half of the width of dimension I of the region. */
/*   WTLENG Integer. */
/*          The number of weights in the basic integration rule. */
/*   G      Real array of dimension (NDIM,WTLENG). */
/*          The fully symmetric sum generators for the rules. */
/*          G(1,J), ..., G(NDIM,J) are the are the generators for the */
/*          points associated with the Jth weights. */
/*   W      Real array of dimension (5,WTLENG). */
/*          The weights for the basic and null rules. */
/*          W(1,1),...,W(1,WTLENG) are weights for the basic rule. */
/*          W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights. */
/*   ERRCOF Real array of dimension 6. */
/*          The error coefficients for the rules. */
/*          It is assumed that the error is computed using: */
/*           IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3) */
/*             THEN ERROR = ERRCOF(3)*N1 */
/*             ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3) */
/*           ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6)) */
/*          where N1-N4 are the null rules, EP is the error */
/*          for the parent */
/*          subregion and ES is the error for the sibling subregion. */
/*   NUMFUN Integer. */
/*          Number of components for the vector integrand. */
/*   FUNSUB Externally declared subroutine. */
/*          For computing the components of the integrand at a point X. */
/*          It must have parameters (NDIM,X,NUMFUN,FUNVLS). */
/*           Input Parameters: */
/*            X      Real array of dimension NDIM. */
/*                   Defines the evaluation point. */
/*            NDIM   Integer. */
/*                   Number of variables for the integrand. */
/*            NUMFUN Integer. */
/*                   Number of components for the vector integrand. */
/*           Output Parameters: */
/*            FUNVLS Real array of dimension NUMFUN. */
/*                   The components of the integrand at the point X. */
/*   SCALES Real array of dimension (3,WTLENG). */
/*          Scaling factors used to construct new null rules based */
/*          on a linear combination of two successive null rules */
/*          in the sequence of null rules. */
/*   NORMS  Real array of dimension (3,WTLENG). */
/*          2**NDIM/(1-norm of the null rule constructed by each of the */
/*          scaling factors.) */
/*   X      Real Array of dimension NDIM. */
/*          A work array. */
/*   NULL   Real array of dimension (NUMFUN, 8) */
/*          A work array. */

/*   ON RETURN */

/*   BASVAL Real array of dimension NUMFUN. */
/*          The values for the basic rule for each component */
/*          of the integrand. */
/*   RGNERR Real array of dimension NUMFUN. */
/*          The error estimates for each component of the integrand. */
/*   DIRECT Real. */
/*          The coordinate axis where the fourth difference of the */
/*          integrand values is largest. */

/* ***REFERENCES */
/*   A.C.Genz and A.A.Malik, An adaptive algorithm for numerical */
/*   integration over an N-dimensional rectangular region, */
/*   J.Comp.Appl.Math., 6:295-302, 1980. */

/*   T.O.Espelid, Integration Rules, Null Rules and Error */
/*   Estimation, Reports in Informatics 33, Dept. of Informatics, */
/*   Univ. of Bergen, 1988. */

/* ***ROUTINES CALLED: DFSHRE, FUNSUB */

/* ***END PROLOGUE DRLHRE */

/*   Global variables. */


/*   Local variables. */


/* ***FIRST EXECUTABLE STATEMENT DRLHRE */


/*       Compute volume of subregion, initialize DIVAXN and rule sums; */
/*       compute fourth differences and new DIVAXN (RGNERR is used */
/*       for a work array here). The integrand values used for the */
/*       fourth divided differences are accumulated in rule arrays. */

    /* Parameter adjustments */
    --x;
    --hwidth;
    --center;
    norms -= 4;
    scales -= 4;
    w -= 6;
    g_dim1 = *ndim;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    --errcof;
    --rgnerr;
    --basval;
    null_dim1 = *numfun;
    null_offset = 1 + null_dim1;
    null -= null_offset;

    /* Function Body */
    rgnvol = 1.;
    divaxn = 1;
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rgnvol *= hwidth[i__];
	x[i__] = center[i__];
	if (hwidth[i__] > hwidth[divaxn]) {
	    divaxn = i__;
	}
/* L10: */
    }
    (*funsub)(ndim, &x[1], numfun, &rgnerr[1]);
    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	basval[j] = w[6] * rgnerr[j];
	for (k = 1; k <= 4; ++k) {
	    null[j + k * null_dim1] = w[k + 6] * rgnerr[j];
/* L20: */
	}
/* L30: */
    }
    difmax = 0.;
/* Computing 2nd power */
    d__1 = g[g_dim1 * 3 + 1] / g[(g_dim1 << 1) + 1];
    ratio = d__1 * d__1;
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = center[i__] - hwidth[i__] * g[(g_dim1 << 1) + 1];
	(*funsub)(ndim, &x[1], numfun, &null[null_dim1 * 5 + 1]);
	x[i__] = center[i__] + hwidth[i__] * g[(g_dim1 << 1) + 1];
	(*funsub)(ndim, &x[1], numfun, &null[null_dim1 * 6 + 1]);
	x[i__] = center[i__] - hwidth[i__] * g[g_dim1 * 3 + 1];
	(*funsub)(ndim, &x[1], numfun, &null[null_dim1 * 7 + 1]);
	x[i__] = center[i__] + hwidth[i__] * g[g_dim1 * 3 + 1];
	(*funsub)(ndim, &x[1], numfun, &null[(null_dim1 << 3) + 1]);
	x[i__] = center[i__];
	difsum = 0.;
	i__2 = *numfun;
	for (j = 1; j <= i__2; ++j) {
	    frthdf = (1 - ratio) * 2 * rgnerr[j] - (null[j + null_dim1 * 7] + 
		    null[j + (null_dim1 << 3)]) + ratio * (null[j + null_dim1 
		    * 5] + null[j + null_dim1 * 6]);

/*           Ignore differences below roundoff */

	    if (rgnerr[j] + frthdf / 4 != rgnerr[j]) {
		difsum += abs(frthdf);
	    }
	    for (k = 1; k <= 4; ++k) {
		null[j + k * null_dim1] = null[j + k * null_dim1] + w[k + 11] 
			* (null[j + null_dim1 * 5] + null[j + null_dim1 * 6]) 
			+ w[k + 16] * (null[j + null_dim1 * 7] + null[j + (
			null_dim1 << 3)]);
/* L40: */
	    }
	    basval[j] = basval[j] + w[11] * (null[j + null_dim1 * 5] + null[j 
		    + null_dim1 * 6]) + w[16] * (null[j + null_dim1 * 7] + 
		    null[j + (null_dim1 << 3)]);
/* L50: */
	}
	if (difsum > difmax) {
	    difmax = difsum;
	    divaxn = i__;
	}
/* L60: */
    }
    *direct = (doublereal) divaxn;

/*    Finish computing the rule values. */

    i__1 = *wtleng;
    for (i__ = 4; i__ <= i__1; ++i__) {
	dfshre_(ndim, &center[1], &hwidth[1], &x[1], &g[i__ * g_dim1 + 1], 
		numfun, (S_fp)funsub, &rgnerr[1], &null[null_dim1 * 5 + 1]);
	i__2 = *numfun;
	for (j = 1; j <= i__2; ++j) {
	    basval[j] += w[i__ * 5 + 1] * rgnerr[j];
	    for (k = 1; k <= 4; ++k) {
		null[j + k * null_dim1] += w[k + 1 + i__ * 5] * rgnerr[j];
/* L70: */
	    }
/* L80: */
	}
/* L90: */
    }

/*    Compute errors. */

    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {

/*    We search for the null rule, in the linear space spanned by two */
/*    successive null rules in our sequence, which gives the greatest */
/*    error estimate among all normalized (1-norm) null rules in this */
/*    space. */

	for (i__ = 1; i__ <= 3; ++i__) {
	    search = 0.;
	    i__2 = *wtleng;
	    for (k = 1; k <= i__2; ++k) {
/* Computing MAX */
		d__2 = search, d__3 = (d__1 = null[j + (i__ + 1) * null_dim1] 
			+ scales[i__ + k * 3] * null[j + i__ * null_dim1], 
			abs(d__1)) * norms[i__ + k * 3];
		search = max(d__2,d__3);
/* L100: */
	    }
	    null[j + i__ * null_dim1] = search;
/* L110: */
	}
	if (errcof[1] * null[j + null_dim1] <= null[j + (null_dim1 << 1)] && 
		errcof[2] * null[j + (null_dim1 << 1)] <= null[j + null_dim1 *
		 3]) {
	    rgnerr[j] = errcof[3] * null[j + null_dim1];
	} else {
/* Computing MAX */
	    d__1 = null[j + null_dim1], d__2 = null[j + (null_dim1 << 1)], 
		    d__1 = max(d__1,d__2), d__2 = null[j + null_dim1 * 3];
	    rgnerr[j] = errcof[4] * max(d__1,d__2);
	}
	rgnerr[j] = rgnvol * rgnerr[j];
	basval[j] = rgnvol * basval[j];
/* L130: */
    }

/* ***END DRLHRE */

    return 0;
} /* drlhre_ */

#ifdef __cplusplus
	}
#endif
