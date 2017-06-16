/* dchhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;

/* Subroutine */ int dchhre_(integer *maxdim, integer *ndim, integer *numfun, 
	integer *mdiv, doublereal *a, doublereal *b, integer *minpts, integer 
	*maxpts, doublereal *epsabs, doublereal *epsrel, integer *key, 
	integer *nw, integer *restar, integer *num, integer *maxsub, integer *
	minsub, integer *keyf, integer *ifail, integer *wtleng)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    //integer pow_ii(integer *, integer *);

    /* Local variables */
    static integer j, limit;

/* ***BEGIN PROLOGUE DCHHRE */
/* ***PURPOSE  DCHHRE checks the validity of the */
/*            input parameters to DCUHRE. */
/* ***DESCRIPTION */
/*            DCHHRE computes NUM, MAXSUB, MINSUB, KEYF, WTLENG and */
/*            IFAIL as functions of the input parameters to DCUHRE, */
/*            and checks the validity of the input parameters to DCUHRE. */

/*   ON ENTRY */

/*     MAXDIM Integer. */
/*            The maximum allowed number of dimensions. */
/*     NDIM   Integer. */
/*            Number of variables. 1 < NDIM <= MAXDIM. */
/*     NUMFUN Integer. */
/*            Number of components of the integral. */
/*     MDIV   Integer. */
/*            MDIV is the number of subregions that are divided in */
/*            each subdivision step in DADHRE. */
/*            MDIV is chosen default to 1. */
/*            For efficient execution on parallel computers */
/*            with NPROC processors MDIV should be set equal to */
/*            the smallest integer such that MOD(2*MDIV,NPROC) = 0. */
/*     A      Real array of dimension NDIM. */
/*            Lower limits of integration. */
/*     B      Real array of dimension NDIM. */
/*            Upper limits of integration. */
/*     MINPTS Integer. */
/*            Minimum number of function evaluations. */
/*     MAXPTS Integer. */
/*            Maximum number of function evaluations. */
/*            The number of function evaluations over each subregion */
/*            is NUM. */
/*            If (KEY = 0 or KEY = 1) and (NDIM = 2) Then */
/*              NUM = 65 */
/*            Elseif (KEY = 0 or KEY = 2) and (NDIM = 3) Then */
/*              NUM = 127 */
/*            Elseif (KEY = 0 and NDIM > 3) or (KEY = 3) Then */
/*              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) + */
/*                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM */
/*            Elseif (KEY = 4) Then */
/*              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM */
/*            MAXPTS >= 3*NUM and MAXPTS >= MINPTS */
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
/*     NW     Integer. */
/*            Defines the length of the working array WORK. */
/*            Let MAXSUB denote the maximum allowed number of subregions */
/*            for the given values of MAXPTS, KEY and NDIM. */
/*            MAXSUB = (MAXPTS-NUM)/(2*NUM) + 1 */
/*            NW should be greater or equal to */
/*            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN + 1 */
/*            For efficient execution on parallel computers */
/*            NW should be chosen greater or equal to */
/*            MAXSUB*(2*NDIM+2*NUMFUN+2) + 17*NUMFUN*MDIV + 1 */
/*            where MDIV is the number of subregions that are divided in */
/*            each subdivision step. */
/*            MDIV is default set internally in DCUHRE equal to 1. */
/*            For efficient execution on parallel computers */
/*            with NPROC processors MDIV should be set equal to */
/*            the smallest integer such that MOD(2*MDIV,NPROC) = 0. */
/*     RESTAR Integer. */
/*            If RESTAR = 0, this is the first attempt to compute */
/*            the integral. */
/*            If RESTAR = 1, then we restart a previous attempt. */

/*   ON RETURN */

/*     NUM    Integer. */
/*            The number of function evaluations over each subregion. */
/*     MAXSUB Integer. */
/*            The maximum allowed number of subregions for the */
/*            given values of MAXPTS, KEY and NDIM. */
/*     MINSUB Integer. */
/*            The minimum allowed number of subregions for the given */
/*            values of MINPTS, KEY and NDIM. */
/*     IFAIL  Integer. */
/*            IFAIL = 0 for normal exit. */
/*            IFAIL = 2 if KEY is less than 0 or KEY greater than 4. */
/*            IFAIL = 3 if NDIM is less than 2 or NDIM greater than */
/*                      MAXDIM. */
/*            IFAIL = 4 if KEY = 1 and NDIM not equal to 2. */
/*            IFAIL = 5 if KEY = 2 and NDIM not equal to 3. */
/*            IFAIL = 6 if NUMFUN less than 1. */
/*            IFAIL = 7 if volume of region of integration is zero. */
/*            IFAIL = 8 if MAXPTS is less than 3*NUM. */
/*            IFAIL = 9 if MAXPTS is less than MINPTS. */
/*            IFAIL = 10 if EPSABS < 0 and EPSREL < 0. */
/*            IFAIL = 11 if NW is too small. */
/*            IFAIL = 12 if unlegal RESTAR. */
/*     KEYF   Integer. */
/*            Key to selected integration rule. */
/*     WTLENG Integer. */
/*            The number of generators of the chosen integration rule. */

/* ***ROUTINES CALLED-NONE */
/* ***END PROLOGUE DCHHRE */

/*   Global variables. */


/*   Local variables. */


/* ***FIRST EXECUTABLE STATEMENT DCHHRE */

    /* Parameter adjustments */
    --b;
    --a;

    /* Function Body */
    *ifail = 0;

/*   Check on legal KEY. */

    if (*key < 0 || *key > 4) {
	*ifail = 2;
	goto L999;
    }

/*   Check on legal NDIM. */

    if (*ndim < 2 || *ndim > *maxdim) {
	*ifail = 3;
	goto L999;
    }

/*   For KEY = 1, NDIM must be equal to 2. */

    if (*key == 1 && *ndim != 2) {
	*ifail = 4;
	goto L999;
    }

/*   For KEY = 2, NDIM must be equal to 3. */

    if (*key == 2 && *ndim != 3) {
	*ifail = 5;
	goto L999;
    }

/*   For KEY = 0, we point at the selected integration rule. */

    if (*key == 0) {
	if (*ndim == 2) {
	    *keyf = 1;
	} else if (*ndim == 3) {
	    *keyf = 2;
	} else {
	    *keyf = 3;
	}
    } else {
	*keyf = *key;
    }

/*   Compute NUM and WTLENG as a function of KEYF and NDIM. */

    if (*keyf == 1) {
	*num = 65;
	*wtleng = 14;
    } else if (*keyf == 2) {
	*num = 127;
	*wtleng = 13;
    } else if (*keyf == 3) {
	*num = (*ndim << 3) + 1 + (*ndim << 1) * (*ndim - 1) + (*ndim << 2) * 
		(*ndim - 1) + (*ndim << 2) * (*ndim - 1) * (*ndim - 2) / 3 + 
		pow_ii(c__2, *ndim);
//		pow_ii(&c__2, ndim);
	*wtleng = 9;
	if (*ndim == 2) {
	    *wtleng = 8;
	}
    } else if (*keyf == 4) {
	*num = *ndim * 6 + 1 + (*ndim << 1) * (*ndim - 1) + pow_ii(c__2, *ndim);
//	*num = *ndim * 6 + 1 + (*ndim << 1) * (*ndim - 1) + pow_ii(&c__2, ndim);
	*wtleng = 6;
    }

/*   Compute MAXSUB. */

    *maxsub = (*maxpts - *num) / (*num << 1) + 1;

/*   Compute MINSUB. */

    *minsub = (*minpts - *num) / (*num << 1) + 1;
    if ((*minpts - *num) % (*num << 1) != 0) {
	++(*minsub);
    }
    *minsub = max(2,*minsub);

/*   Check on positive NUMFUN. */

    if (*numfun < 1) {
	*ifail = 6;
	goto L999;
    }

/*   Check on legal upper and lower limits of integration. */

    i__1 = *ndim;
    for (j = 1; j <= i__1; ++j) {
	if (a[j] - b[j] == 0.) {
	    *ifail = 7;
	    goto L999;
	}
/* L10: */
    }

/*   Check on MAXPTS < 3*NUM. */

    if (*maxpts < *num * 3) {
	*ifail = 8;
	goto L999;
    }

/*   Check on MAXPTS >= MINPTS. */

    if (*maxpts < *minpts) {
	*ifail = 9;
	goto L999;
    }

/*   Check on legal accuracy requests. */

    if (*epsabs < 0. && *epsrel < 0.) {
	*ifail = 10;
	goto L999;
    }

/*   Check on big enough double precision workspace. */

    limit = *maxsub * ((*ndim << 1) + (*numfun << 1) + 2) + *mdiv * 17 * *
	    numfun + 1;
    if (*nw < limit) {
	*ifail = 11;
	goto L999;
    }

/*    Check on legal RESTAR. */

    if (*restar != 0 && *restar != 1) {
	*ifail = 12;
	goto L999;
    }
L999:
    return 0;

/* ***END DCHHRE */

} /* dchhre_ */

#ifdef __cplusplus
	}
#endif
