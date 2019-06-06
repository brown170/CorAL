/* dinhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int dinhre_(integer *ndim, integer *key, integer *wtleng,
	doublereal *w, doublereal *g, doublereal *errcof, doublereal *rulpts,
	doublereal *scales, doublereal *norms)
{
    /* System generated locals */
    integer g_dim1, g_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    //integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer i__, j, k;
    static doublereal we[14];
    extern /* Subroutine */ int d113re_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *), d132re_(integer *, doublereal *,
	    doublereal *, doublereal *, doublereal *), d07hre_(integer *,
	    integer *, doublereal *, doublereal *, doublereal *, doublereal *)
	    , d09hre_(integer *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *);

/* ***BEGIN PROLOGUE DINHRE */
/* ***PURPOSE DINHRE computes abscissas and weights of the integration */
/*            rule and the null rules to be used in error estimation. */
/*            These are computed as functions of NDIM and KEY. */
/* ***DESCRIPTION DINHRE will for given values of NDIM and KEY compute or */
/*            select the correct values of the abscissas and */
/*            corresponding weights for different */
/*            integration rules and null rules and assign them to */
/*            G and W. */
/*            The heuristic error coefficients ERRCOF */
/*            will be computed as a function of KEY. */
/*            Scaling factors SCALES and normalization factors NORMS */
/*            used in the error estimation are computed. */


/*   ON ENTRY */

/*     NDIM   Integer. */
/*            Number of variables. */
/*     KEY    Integer. */
/*            Key to selected local integration rule. */
/*     WTLENG Integer. */
/*            The number of weights in each of the rules. */

/*   ON RETURN */

/*     W      Real array of dimension (5,WTLENG). */
/*            The weights for the basic and null rules. */
/*            W(1,1), ...,W(1,WTLENG) are weights for the basic rule. */
/*            W(I,1), ...,W(I,WTLENG), for I > 1 are null rule weights. */
/*     G      Real array of dimension (NDIM,WTLENG). */
/*            The fully symmetric sum generators for the rules. */
/*            G(1,J),...,G(NDIM,J) are the generators for the points */
/*            associated with the the Jth weights. */
/*     ERRCOF Real array of dimension 6. */
/*            Heuristic error coefficients that are used in the */
/*            error estimation in BASRUL. */
/*            It is assumed that the error is computed using: */
/*             IF (N1*ERRCOF(1) < N2 and N2*ERRCOF(2) < N3) */
/*               THEN ERROR = ERRCOF(3)*N1 */
/*               ELSE ERROR = ERRCOF(4)*MAX(N1, N2, N3) */
/*             ERROR = ERROR + EP*(ERRCOF(5)*ERROR/(ES+ERROR)+ERRCOF(6)) */
/*            where N1-N3 are the null rules, EP is the error for */
/*            the parent */
/*            subregion and ES is the error for the sibling subregion. */
/*     RULPTS Real array of dimension WTLENG. */
/*            A work array containing the number of points produced by */
/*            each generator of the selected rule. */
/*     SCALES Real array of dimension (3,WTLENG). */
/*            Scaling factors used to construct new null rules, */
/*            N1, N2 and N3, */
/*            based on a linear combination of two successive null rules */
/*            in the sequence of null rules. */
/*     NORMS  Real array of dimension (3,WTLENG). */
/*            2**NDIM/(1-norm of the null rule constructed by each of */
/*            the scaling factors.) */

/* ***ROUTINES CALLED  D132RE,D113RE,D07HRE,D09HRE */
/* ***END PROLOGUE DINHRE */

/*   Global variables. */


/*   Local variables. */


/* ***FIRST EXECUTABLE STATEMENT DINHRE */

/*   Compute W, G and ERRCOF. */

    /* Parameter adjustments */
    norms -= 4;
    scales -= 4;
    --rulpts;
    g_dim1 = *ndim;
    g_offset = 1 + g_dim1;
    g -= g_offset;
    w -= 6;
    --errcof;

    /* Function Body */
    if (*key == 1) {
	d132re_(wtleng, &w[6], &g[g_offset], &errcof[1], &rulpts[1]);
    } else if (*key == 2) {
	d113re_(wtleng, &w[6], &g[g_offset], &errcof[1], &rulpts[1]);
    } else if (*key == 3) {
	d09hre_(ndim, wtleng, &w[6], &g[g_offset], &errcof[1], &rulpts[1]);
    } else if (*key == 4) {
	d07hre_(ndim, wtleng, &w[6], &g[g_offset], &errcof[1], &rulpts[1]);
    }

/*   Compute SCALES and NORMS. */

    for (k = 1; k <= 3; ++k) {
	i__1 = *wtleng;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (w[k + 1 + i__ * 5] != 0.) {
		scales[k + i__ * 3] = -w[k + 2 + i__ * 5] / w[k + 1 + i__ * 5]
			;
	    } else {
		scales[k + i__ * 3] = 100.;
	    }
	    i__2 = *wtleng;
	    for (j = 1; j <= i__2; ++j) {
		we[j - 1] = w[k + 2 + j * 5] + scales[k + i__ * 3] * w[k + 1
			+ j * 5];
/* L30: */
	    }
	    norms[k + i__ * 3] = 0.;
	    i__2 = *wtleng;
	    for (j = 1; j <= i__2; ++j) {
		norms[k + i__ * 3] += rulpts[j] * (d__1 = we[j - 1], ABS(d__1)
			);
/* L40: */
	    }
	    norms[k + i__ * 3] = pow_ii(c__2, *ndim) / norms[k + i__ * 3];
	    //norms[k + i__ * 3] = pow_ii(&c__2, ndim) / norms[k + i__ * 3];
/* L50: */
	}
/* L100: */
    }
    return 0;

/* ***END DINHRE */

} /* dinhre_ */

#ifdef __cplusplus
	}
#endif
