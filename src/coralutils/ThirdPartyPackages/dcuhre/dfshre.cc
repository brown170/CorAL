/* dfshre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dfshre_(integer *ndim, doublereal *center, doublereal *
	hwidth, doublereal *x, doublereal *g, integer *numfun, S_fp funsub, 
	doublereal *fulsms, doublereal *funvls)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, j, l;
    static doublereal gi, gl;
    static integer ixchng, lxchng;

/* ***BEGIN PROLOGUE DFSHRE */
/* ***KEYWORDS fully symmetric sum */
/* ***PURPOSE  To compute fully symmetric basic rule sums */
/* ***AUTHOR   Alan Genz, Computer Science Department, Washington */
/*            State University, Pullman, WA 99163-1210 USA */
/* ***LAST MODIFICATION 88-04-08 */
/* ***DESCRIPTION DFSHRE computes a fully symmetric sum for a vector */
/*            of integrand values over a hyper-rectangular region. */
/*            The sum is fully symmetric with respect to the center of */
/*            the region and is taken over all sign changes and */
/*            permutations of the generators for the sum. */

/*   ON ENTRY */

/*   NDIM   Integer. */
/*          Number of variables. */
/*   CENTER Real array of dimension NDIM. */
/*          The coordinates for the center of the region. */
/*   HWIDTH Real Array of dimension NDIM. */
/*          HWIDTH(I) is half of the width of dimension I of the region. */
/*   X      Real Array of dimension NDIM. */
/*          A work array. */
/*   G      Real Array of dimension NDIM. */
/*          The generators for the fully symmetric sum. These MUST BE */
/*          non-negative and non-increasing. */
/*   NUMFUN Integer. */
/*          Number of components for the vector integrand. */
/*   FUNSUB Externally declared subroutine. */
/*          For computing the components of the integrand at a point X. */
/*          It must have parameters (NDIM, X, NUMFUN, FUNVLS). */
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
/*   ON RETURN */

/*   FULSMS Real array of dimension NUMFUN. */
/*          The values for the fully symmetric sums for each component */
/*          of the integrand. */
/*   FUNVLS Real array of dimension NUMFUN. */
/*          A work array. */

/* ***ROUTINES CALLED: FUNSUB */

/* ***END PROLOGUE DFSHRE */

/*   Global variables. */


/*   Local variables. */


/* ***FIRST EXECUTABLE STATEMENT DFSHRE */

    /* Parameter adjustments */
    --g;
    --x;
    --hwidth;
    --center;
    --funvls;
    --fulsms;

    /* Function Body */
    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	fulsms[j] = 0.;
/* L10: */
    }

/*     Compute centrally symmetric sum for permutation of G */

L20:
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = center[i__] + g[i__] * hwidth[i__];
/* L30: */
    }
L40:
    (*funsub)(ndim, &x[1], numfun, &funvls[1]);
    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	fulsms[j] += funvls[j];
/* L50: */
    }
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = -g[i__];
	x[i__] = center[i__] + g[i__] * hwidth[i__];
	if (g[i__] < 0.) {
	    goto L40;
	}
/* L60: */
    }

/*       Find next distinct permutation of G and loop back for next sum. */
/*       Permutations are generated in reverse lexicographic order. */

    i__1 = *ndim;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (g[i__ - 1] > g[i__]) {
	    gi = g[i__];
	    ixchng = i__ - 1;
	    i__2 = (i__ - 1) / 2;
	    for (l = 1; l <= i__2; ++l) {
		gl = g[l];
		g[l] = g[i__ - l];
		g[i__ - l] = gl;
		if (gl <= gi) {
		    --ixchng;
		}
		if (g[l] > gi) {
		    lxchng = l;
		}
/* L70: */
	    }
	    if (g[ixchng] <= gi) {
		ixchng = lxchng;
	    }
	    g[i__] = g[ixchng];
	    g[ixchng] = gi;
	    goto L20;
	}
/* L80: */
    }

/*     Restore original order to generators */

    i__1 = *ndim / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = g[i__];
	g[i__] = g[*ndim - i__ + 1];
	g[*ndim - i__ + 1] = gi;
/* L90: */
    }

/* ***END DFSHRE */

    return 0;
} /* dfshre_ */

#ifdef __cplusplus
	}
#endif
