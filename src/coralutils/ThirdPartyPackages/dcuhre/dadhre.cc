/* dadhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* Subroutine */ int dadhre_(integer *ndim, integer *numfun, integer *mdiv,
	doublereal *a, doublereal *b, integer *minsub, integer *maxsub, U_fp
	funsub, doublereal *epsabs, doublereal *epsrel, integer *key, integer
	*restar, integer *num, integer *lenw, integer *wtleng, doublereal *
	result, doublereal *abserr, integer *neval, integer *nsub, integer *
	ifail, doublereal *values, doublereal *errors, doublereal *centrs,
	doublereal *hwidts, doublereal *greate, doublereal *dir, doublereal *
	oldres, doublereal *work, doublereal *g, doublereal *w, doublereal *
	rulpts, doublereal *center, doublereal *hwidth, doublereal *x,
	doublereal *scales, doublereal *norms)
{
    /* System generated locals */
    integer values_dim1, values_offset, errors_dim1, errors_offset,
	    centrs_dim1, centrs_offset, hwidts_dim1, hwidts_offset,
	    oldres_dim1, oldres_offset, g_dim1, g_dim2, g_offset, x_dim1,
	    x_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    static integer i__, j, k, l1;
    static doublereal est1, est2;
    static integer ndiv, index;
    static doublereal oldcen;
    extern /* Subroutine */ int dinhre_(integer *, integer *, integer *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *);
    static integer direct;
    static doublereal errcof[6];
    extern /* Subroutine */ int drlhre_(integer *, doublereal *, doublereal *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *,
	    U_fp, doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *), dtrhre_(integer *,
	    integer *, integer *, integer *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *,
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer intsgn, sbrgns, pointr;

/* ***BEGIN PROLOGUE DADHRE */
/* ***KEYWORDS automatic multidimensional integrator, */
/*            n-dimensional hyper-rectangles, */
/*            general purpose, global adaptive */
/* ***PURPOSE  The routine calculates an approximation to a given */
/*            vector of definite integrals, I, over a hyper-rectangular */
/*            region hopefully satisfying for each component of I the */
/*            following claim for accuracy: */
/*            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K))) */
/* ***DESCRIPTION Computation of integrals over hyper-rectangular */
/*            regions. */
/*            DADHRE repeatedly subdivides the region */
/*            of integration and estimates the integrals and the */
/*            errors over the subregions with  greatest */
/*            estimated errors until the error request */
/*            is met or MAXSUB subregions are stored. */
/*            The regions are divided in two equally sized parts along */
/*            the direction with greatest absolute fourth divided */
/*            difference. */

/*   ON ENTRY */

/*     NDIM   Integer. */
/*            Number of variables. 1 < NDIM <= MAXDIM. */
/*     NUMFUN Integer. */
/*            Number of components of the integral. */
/*     MDIV   Integer. */
/*            Defines the number of new subregions that are divided */
/*            in each subdivision step. */
/*     A      Real array of dimension NDIM. */
/*            Lower limits of integration. */
/*     B      Real array of dimension NDIM. */
/*            Upper limits of integration. */
/*     MINSUB Integer. */
/*            The computations proceed until there are at least */
/*            MINSUB subregions in the data structure. */
/*     MAXSUB Integer. */
/*            The computations proceed until there are at most */
/*            MAXSUB subregions in the data structure. */

/*     FUNSUB Externally declared subroutine for computing */
/*            all components of the integrand in the given */
/*            evaluation point. */
/*            It must have parameters (NDIM,X,NUMFUN,FUNVLS) */
/*            Input parameters: */
/*              NDIM   Integer that defines the dimension of the */
/*                     integral. */
/*              X      Real array of dimension NDIM */
/*                     that defines the evaluation point. */
/*              NUMFUN Integer that defines the number of */
/*                     components of I. */
/*            Output parameter: */
/*              FUNVLS Real array of dimension NUMFUN */
/*                     that defines NUMFUN components of the integrand. */

/*     EPSABS Real. */
/*            Requested absolute error. */
/*     EPSREL Real. */
/*            Requested relative error. */
/*     KEY    Integer. */
/*            Key to selected local integration rule. */
/*            KEY = 0 is the default value. */
/*                  For NDIM = 2 the degree 13 rule is selected. */
/*                  For NDIM = 3 the degree 11 rule is selected. */
/*                  For NDIM > 3 the degree  9 rule is selected. */
/*            KEY = 1 gives the user the 2 dimensional degree 13 */
/*                  integration rule that uses 65 evaluation points. */
/*            KEY = 2 gives the user the 3 dimensional degree 11 */
/*                  integration rule that uses 127 evaluation points. */
/*            KEY = 3 gives the user the degree 9 integration rule. */
/*            KEY = 4 gives the user the degree 7 integration rule. */
/*                  This is the recommended rule for problems that */
/*                  require great adaptivity. */
/*     RESTAR Integer. */
/*            If RESTAR = 0, this is the first attempt to compute */
/*            the integral. */
/*            If RESTAR = 1, then we restart a previous attempt. */
/*            (In this case the output parameters from DADHRE */
/*            must not be changed since the last */
/*            exit from DADHRE.) */
/*     NUM    Integer. */
/*            The number of function evaluations over each subregion. */
/*     LENW   Integer. */
/*            Defines the length of the working array WORK. */
/*            LENW should be greater or equal to */
/*            16*MDIV*NUMFUN. */
/*     WTLENG Integer. */
/*            The number of weights in the basic integration rule. */
/*     NSUB   Integer. */
/*            If RESTAR = 1, then NSUB must specify the number */
/*            of subregions stored in the previous call to DADHRE. */

/*   ON RETURN */

/*     RESULT Real array of dimension NUMFUN. */
/*            Approximations to all components of the integral. */
/*     ABSERR Real array of dimension NUMFUN. */
/*            Estimates of absolute accuracies. */
/*     NEVAL  Integer. */
/*            Number of function evaluations used by DADHRE. */
/*     NSUB   Integer. */
/*            Number of stored subregions. */
/*     IFAIL  Integer. */
/*            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or */
/*              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXSUB or less */
/*              subregions processed for all values of K, */
/*              1 <=  K <=  NUMFUN. */
/*            IFAIL = 1 if MAXSUB was too small for DADHRE */
/*              to obtain the required accuracy. In this case DADHRE */
/*              returns values of RESULT with estimated absolute */
/*              accuracies ABSERR. */
/*     VALUES Real array of dimension (NUMFUN,MAXSUB). */
/*            Used to store estimated values of the integrals */
/*            over the subregions. */
/*     ERRORS Real array of dimension (NUMFUN,MAXSUB). */
/*            Used to store the corresponding estimated errors. */
/*     CENTRS Real array of dimension (NDIM,MAXSUB). */
/*            Used to store the centers of the stored subregions. */
/*     HWIDTS Real array of dimension (NDIM,MAXSUB). */
/*            Used to store the half widths of the stored subregions. */
/*     GREATE Real array of dimension MAXSUB. */
/*            Used to store the greatest estimated errors in */
/*            all subregions. */
/*     DIR    Real array of dimension MAXSUB. */
/*            DIR is used to store the directions for */
/*            further subdivision. */
/*     OLDRES Real array of dimension (NUMFUN,MDIV). */
/*            Used to store old estimates of the integrals over the */
/*            subregions. */
/*     WORK   Real array of dimension LENW. */
/*            Used  in DRLHRE and DTRHRE. */
/*     G      Real array of dimension (NDIM,WTLENG,2*MDIV). */
/*            The fully symmetric sum generators for the rules. */
/*            G(1,J,1),...,G(NDIM,J,1) are the generators for the */
/*            points associated with the Jth weights. */
/*            When MDIV subregions are divided in 2*MDIV */
/*            subregions, the subregions may be processed on different */
/*            processors and we must make a copy of the generators */
/*            for each processor. */
/*     W      Real array of dimension (5,WTLENG). */
/*            The weights for the basic and null rules. */
/*            W(1,1), ..., W(1,WTLENG) are weights for the basic rule. */
/*            W(I,1), ..., W(I,WTLENG) , for I > 1 are null rule weights. */
/*     RULPTS Real array of dimension WTLENG. */
/*            Work array used in DINHRE. */
/*     CENTER Real array of dimension NDIM. */
/*            Work array used in DTRHRE. */
/*     HWIDTH Real array of dimension NDIM. */
/*            Work array used in DTRHRE. */
/*     X      Real array of dimension (NDIM,2*MDIV). */
/*            Work array used in DRLHRE. */
/*     SCALES Real array of dimension (3,WTLENG). */
/*            Work array used by DINHRE and DRLHRE. */
/*     NORMS  Real array of dimension (3,WTLENG). */
/*            Work array used by DINHRE and DRLHRE. */

/* ***REFERENCES */

/*   P. van Dooren and L. de Ridder, Algorithm 6, An adaptive algorithm */
/*   for numerical integration over an n-dimensional cube, J.Comput.Appl. */
/*   Math. 2(1976)207-217. */

/*   A.C.Genz and A.A.Malik, Algorithm 019. Remarks on algorithm 006: */
/*   An adaptive algorithm for numerical integration over an */
/*   N-dimensional rectangular region,J.Comput.Appl.Math. 6(1980)295-302. */

/* ***ROUTINES CALLED DTRHRE,DINHRE,DRLHRE */
/* ***END PROLOGUE DADHRE */

/*   Global variables. */


/*   Local variables. */

/*   INTSGN is used to get correct sign on the integral. */
/*   SBRGNS is the number of stored subregions. */
/*   NDIV   The number of subregions to be divided in each main step. */
/*   POINTR Pointer to the position in the data structure where */
/*          the new subregions are to be stored. */
/*   DIRECT Direction of subdivision. */
/*   ERRCOF Heuristic error coeff. defined in DINHRE and used by DRLHRE */
/*          and DADHRE. */


/* ***FIRST EXECUTABLE STATEMENT DADHRE */

/*   Get the correct sign on the integral. */

    /* Parameter adjustments */
    --hwidth;
    --center;
    --b;
    --a;
    --abserr;
    --result;
    x_dim1 = *ndim;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    oldres_dim1 = *numfun;
    oldres_offset = 1 + oldres_dim1;
    oldres -= oldres_offset;
    --dir;
    --greate;
    hwidts_dim1 = *ndim;
    hwidts_offset = 1 + hwidts_dim1;
    hwidts -= hwidts_offset;
    centrs_dim1 = *ndim;
    centrs_offset = 1 + centrs_dim1;
    centrs -= centrs_offset;
    errors_dim1 = *numfun;
    errors_offset = 1 + errors_dim1;
    errors -= errors_offset;
    values_dim1 = *numfun;
    values_offset = 1 + values_dim1;
    values -= values_offset;
    --work;
    norms -= 4;
    scales -= 4;
    --rulpts;
    w -= 6;
    g_dim1 = *ndim;
    g_dim2 = *wtleng;
    g_offset = 1 + g_dim1 * (1 + g_dim2);
    g -= g_offset;

    /* Function Body */
    intsgn = 1;
    i__1 = *ndim;
    for (j = 1; j <= i__1; ++j) {
	if (b[j] < a[j]) {
	    intsgn = -intsgn;
	}
/* L10: */
    }

/*   Call DINHRE to compute the weights and abscissas of */
/*   the function evaluation points. */

    dinhre_(ndim, key, wtleng, &w[6], &g[g_offset], errcof, &rulpts[1], &
	    scales[4], &norms[4]);

/*   If RESTAR = 1, then this is a restart run. */

    if (*restar == 1) {
	sbrgns = *nsub;
	goto L110;
    }

/*   Initialize the SBRGNS, CENTRS and HWIDTS. */

    sbrgns = 1;
    i__1 = *ndim;
    for (j = 1; j <= i__1; ++j) {
	centrs[j + centrs_dim1] = (a[j] + b[j]) / 2;
	hwidts[j + hwidts_dim1] = (d__1 = b[j] - a[j], ABS(d__1)) / 2;
/* L15: */
    }

/*   Initialize RESULT, ABSERR and NEVAL. */

    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	result[j] = 0.;
	abserr[j] = 0.;
/* L20: */
    }
    *neval = 0;

/*   Apply DRLHRE over the whole region. */

    drlhre_(ndim, &centrs[centrs_dim1 + 1], &hwidts[hwidts_dim1 + 1], wtleng,
	    &g[g_offset], &w[6], errcof, numfun, (U_fp)funsub, &scales[4], &
	    norms[4], &x[x_offset], &work[1], &values[values_dim1 + 1], &
	    errors[errors_dim1 + 1], &dir[1]);
    *neval += *num;

/*   Add the computed values to RESULT and ABSERR. */

    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	result[j] += values[j + values_dim1];
/* L55: */
    }
    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	abserr[j] += errors[j + errors_dim1];
/* L65: */
    }

/*   Store results in heap. */

    index = 1;
    dtrhre_(&c__2, ndim, numfun, &index, &values[values_offset], &errors[
	    errors_offset], &centrs[centrs_offset], &hwidts[hwidts_offset], &
	    greate[1], &work[1], &work[*numfun + 1], &center[1], &hwidth[1], &
	    dir[1]);

/* ***End initialisation. */

/* ***Begin loop while the error is too great, */
/*   and SBRGNS+1 is less than MAXSUB. */

L110:
    if (sbrgns + 1 <= *maxsub) {

/*   If we are allowed to divide further, */
/*   prepare to apply basic rule over each half of the */
/*   NDIV subregions with greatest errors. */
/*   If MAXSUB is great enough, NDIV = MDIV */

	if (*mdiv > 1) {
	    ndiv = *maxsub - sbrgns;
/* Computing MIN */
	    i__1 = min(ndiv,*mdiv);
	    ndiv = min(i__1,sbrgns);
	} else {
	    ndiv = 1;
	}

/*   Divide the NDIV subregions in two halves, and compute */
/*   integral and error over each half. */

	i__1 = ndiv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    pointr = sbrgns + ndiv + 1 - i__;

/*   Adjust RESULT and ABSERR. */

	    i__2 = *numfun;
	    for (j = 1; j <= i__2; ++j) {
		result[j] -= values[j + values_dim1];
		abserr[j] -= errors[j + errors_dim1];
/* L115: */
	    }

/*   Compute first half region. */

	    i__2 = *ndim;
	    for (j = 1; j <= i__2; ++j) {
		centrs[j + pointr * centrs_dim1] = centrs[j + centrs_dim1];
		hwidts[j + pointr * hwidts_dim1] = hwidts[j + hwidts_dim1];
/* L120: */
	    }
	    direct = (integer) dir[1];
	    dir[pointr] = (doublereal) direct;
	    hwidts[direct + pointr * hwidts_dim1] = hwidts[direct +
		    hwidts_dim1] / 2;
	    oldcen = centrs[direct + centrs_dim1];
	    centrs[direct + pointr * centrs_dim1] = oldcen - hwidts[direct +
		    pointr * hwidts_dim1];

/*   Save the computed values of the integrals. */

	    i__2 = *numfun;
	    for (j = 1; j <= i__2; ++j) {
		oldres[j + (ndiv - i__ + 1) * oldres_dim1] = values[j +
			values_dim1];
/* L125: */
	    }

/*   Adjust the heap. */

	    dtrhre_(&c__1, ndim, numfun, &sbrgns, &values[values_offset], &
		    errors[errors_offset], &centrs[centrs_offset], &hwidts[
		    hwidts_offset], &greate[1], &work[1], &work[*numfun + 1],
		    &center[1], &hwidth[1], &dir[1]);

/*   Compute second half region. */

	    i__2 = *ndim;
	    for (j = 1; j <= i__2; ++j) {
		centrs[j + (pointr - 1) * centrs_dim1] = centrs[j + pointr *
			centrs_dim1];
		hwidts[j + (pointr - 1) * hwidts_dim1] = hwidts[j + pointr *
			hwidts_dim1];
/* L130: */
	    }
	    centrs[direct + (pointr - 1) * centrs_dim1] = oldcen + hwidts[
		    direct + pointr * hwidts_dim1];
	    hwidts[direct + (pointr - 1) * hwidts_dim1] = hwidts[direct +
		    pointr * hwidts_dim1];
	    dir[pointr - 1] = (doublereal) direct;
/* L150: */
	}

/*   Make copies of the generators for each processor. */

	i__1 = ndiv << 1;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    i__2 = *ndim;
	    for (j = 1; j <= i__2; ++j) {
		i__3 = *wtleng;
		for (k = 1; k <= i__3; ++k) {
		    g[j + (k + i__ * g_dim2) * g_dim1] = g[j + (k + g_dim2) *
			    g_dim1];
/* L190: */
		}
	    }
	}

/*   Apply basic rule. */

/* vd$l cncall */
	i__3 = ndiv << 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    index = sbrgns + i__;
	    l1 = (i__ - 1 << 3) * *numfun + 1;
	    drlhre_(ndim, &centrs[index * centrs_dim1 + 1], &hwidts[index *
		    hwidts_dim1 + 1], wtleng, &g[(i__ * g_dim2 + 1) * g_dim1
		    + 1], &w[6], errcof, numfun, (U_fp)funsub, &scales[4], &
		    norms[4], &x[i__ * x_dim1 + 1], &work[l1], &values[index *
		     values_dim1 + 1], &errors[index * errors_dim1 + 1], &dir[
		    index]);
/* L200: */
	}
	*neval += (ndiv << 1) * *num;

/*   Add new contributions to RESULT. */

	i__3 = ndiv << 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    i__2 = *numfun;
	    for (j = 1; j <= i__2; ++j) {
		result[j] += values[j + (sbrgns + i__) * values_dim1];
/* L210: */
	    }
/* L220: */
	}

/*   Check consistency of results and if necessary adjust */
/*   the estimated errors. */

	i__3 = ndiv;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    greate[sbrgns + (i__ << 1) - 1] = 0.;
	    greate[sbrgns + (i__ << 1)] = 0.;
	    i__2 = *numfun;
	    for (j = 1; j <= i__2; ++j) {
		est1 = (d__1 = oldres[j + i__ * oldres_dim1] - (values[j + (
			sbrgns + (i__ << 1) - 1) * values_dim1] + values[j + (
			sbrgns + (i__ << 1)) * values_dim1]), ABS(d__1));
		est2 = errors[j + (sbrgns + (i__ << 1) - 1) * errors_dim1] +
			errors[j + (sbrgns + (i__ << 1)) * errors_dim1];
		if (est2 > 0.) {
		    errors[j + (sbrgns + (i__ << 1) - 1) * errors_dim1] *=
			    errcof[4] * est1 / est2 + 1;
		    errors[j + (sbrgns + (i__ << 1)) * errors_dim1] *= errcof[
			    4] * est1 / est2 + 1;
		}
		errors[j + (sbrgns + (i__ << 1) - 1) * errors_dim1] += errcof[
			5] * est1;
		errors[j + (sbrgns + (i__ << 1)) * errors_dim1] += errcof[5] *
			 est1;
		if (errors[j + (sbrgns + (i__ << 1) - 1) * errors_dim1] >
			greate[sbrgns + (i__ << 1) - 1]) {
		    greate[sbrgns + (i__ << 1) - 1] = errors[j + (sbrgns + (
			    i__ << 1) - 1) * errors_dim1];
		}
		if (errors[j + (sbrgns + (i__ << 1)) * errors_dim1] > greate[
			sbrgns + (i__ << 1)]) {
		    greate[sbrgns + (i__ << 1)] = errors[j + (sbrgns + (i__ <<
			     1)) * errors_dim1];
		}
		abserr[j] = abserr[j] + errors[j + (sbrgns + (i__ << 1) - 1) *
			 errors_dim1] + errors[j + (sbrgns + (i__ << 1)) *
			errors_dim1];
/* L230: */
	    }
/* L240: */
	}

/*   Store results in heap. */

	i__3 = ndiv << 1;
	for (i__ = 1; i__ <= i__3; ++i__) {
	    index = sbrgns + i__;
	    dtrhre_(&c__2, ndim, numfun, &index, &values[values_offset], &
		    errors[errors_offset], &centrs[centrs_offset], &hwidts[
		    hwidts_offset], &greate[1], &work[1], &work[*numfun + 1],
		    &center[1], &hwidth[1], &dir[1]);
/* L250: */
	}
	sbrgns += ndiv << 1;

/*   Check for termination. */

	if (sbrgns < *minsub) {
	    goto L110;
	}
	i__3 = *numfun;
	for (j = 1; j <= i__3; ++j) {
	    if (abserr[j] > *epsrel * (d__1 = result[j], ABS(d__1)) && abserr[
		    j] > *epsabs) {
		goto L110;
	    }
/* L255: */
	}
	*ifail = 0;
	goto L499;

/*   Else we did not succeed with the */
/*   given value of MAXSUB. */

    } else {
	*ifail = 1;
    }

/*   Compute more accurate values of RESULT and ABSERR. */

L499:
    i__3 = *numfun;
    for (j = 1; j <= i__3; ++j) {
	result[j] = 0.;
	abserr[j] = 0.;
/* L500: */
    }
    i__3 = sbrgns;
    for (i__ = 1; i__ <= i__3; ++i__) {
	i__2 = *numfun;
	for (j = 1; j <= i__2; ++j) {
	    result[j] += values[j + i__ * values_dim1];
	    abserr[j] += errors[j + i__ * errors_dim1];
/* L505: */
	}
/* L510: */
    }

/*   Compute correct sign on the integral. */

    i__3 = *numfun;
    for (j = 1; j <= i__3; ++j) {
	result[j] *= intsgn;
/* L600: */
    }
    *nsub = sbrgns;
    return 0;

/* ***END DADHRE */

} /* dadhre_ */

#ifdef __cplusplus
	}
#endif
