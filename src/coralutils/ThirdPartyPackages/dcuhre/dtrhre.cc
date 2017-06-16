/* dtrhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int dtrhre_(integer *dvflag, integer *ndim, integer *numfun, 
	integer *sbrgns, doublereal *values, doublereal *errors, doublereal *
	centrs, doublereal *hwidts, doublereal *greate, doublereal *error, 
	doublereal *value, doublereal *center, doublereal *hwidth, doublereal 
	*dir)
{
    /* System generated locals */
    integer values_dim1, values_offset, errors_dim1, errors_offset, 
	    centrs_dim1, centrs_offset, hwidts_dim1, hwidts_offset, i__1;

    /* Local variables */
    static integer j;
    static doublereal great, direct;
    static integer subrgn, subtmp;

/* ***BEGIN PROLOGUE DTRHRE */
/* ***PURPOSE DTRHRE maintains a heap of subregions. */
/* ***DESCRIPTION DTRHRE maintains a heap of subregions. */
/*            The subregions are ordered according to the size */
/*            of the greatest error estimates of each subregion(GREATE). */

/*   PARAMETERS */

/*     DVFLAG Integer. */
/*            If DVFLAG = 1, we remove the subregion with */
/*            greatest error from the heap. */
/*            If DVFLAG = 2, we insert a new subregion in the heap. */
/*     NDIM   Integer. */
/*            Number of variables. */
/*     NUMFUN Integer. */
/*            Number of components of the integral. */
/*     SBRGNS Integer. */
/*            Number of subregions in the heap. */
/*     VALUES Real array of dimension (NUMFUN,SBRGNS). */
/*            Used to store estimated values of the integrals */
/*            over the subregions. */
/*     ERRORS Real array of dimension (NUMFUN,SBRGNS). */
/*            Used to store the corresponding estimated errors. */
/*     CENTRS Real array of dimension (NDIM,SBRGNS). */
/*            Used to store the center limits of the stored subregions. */
/*     HWIDTS Real array of dimension (NDIM,SBRGNS). */
/*            Used to store the hwidth limits of the stored subregions. */
/*     GREATE Real array of dimension SBRGNS. */
/*            Used to store the greatest estimated errors in */
/*            all subregions. */
/*     ERROR  Real array of dimension NUMFUN. */
/*            Used as intermediate storage for the error of a subregion. */
/*     VALUE  Real array of dimension NUMFUN. */
/*            Used as intermediate storage for the estimate */
/*            of the integral over a subregion. */
/*     CENTER Real array of dimension NDIM. */
/*            Used as intermediate storage for the center of */
/*            the subregion. */
/*     HWIDTH Real array of dimension NDIM. */
/*            Used as intermediate storage for the half width of */
/*            the subregion. */
/*     DIR    Integer array of dimension SBRGNS. */
/*            DIR is used to store the directions for */
/*            further subdivision. */

/* ***ROUTINES CALLED-NONE */
/* ***END PROLOGUE DTRHRE */

/*   Global variables. */


/*   Local variables. */

/*   GREAT  is used as intermediate storage for the greatest error of a */
/*          subregion. */
/*   DIRECT is used as intermediate storage for the direction of further */
/*          subdivision. */
/*   SUBRGN Position of child/parent subregion in the heap. */
/*   SUBTMP Position of parent/child subregion in the heap. */


/* ***FIRST EXECUTABLE STATEMENT DTRHRE */

/*   Save values to be stored in their correct place in the heap. */

    /* Parameter adjustments */
    --hwidth;
    --center;
    hwidts_dim1 = *ndim;
    hwidts_offset = 1 + hwidts_dim1;
    hwidts -= hwidts_offset;
    centrs_dim1 = *ndim;
    centrs_offset = 1 + centrs_dim1;
    centrs -= centrs_offset;
    --value;
    --error;
    errors_dim1 = *numfun;
    errors_offset = 1 + errors_dim1;
    errors -= errors_offset;
    values_dim1 = *numfun;
    values_offset = 1 + values_dim1;
    values -= values_offset;
    --greate;
    --dir;

    /* Function Body */
    great = greate[*sbrgns];
    direct = dir[*sbrgns];
    i__1 = *numfun;
    for (j = 1; j <= i__1; ++j) {
	error[j] = errors[j + *sbrgns * errors_dim1];
	value[j] = values[j + *sbrgns * values_dim1];
/* L5: */
    }
    i__1 = *ndim;
    for (j = 1; j <= i__1; ++j) {
	center[j] = centrs[j + *sbrgns * centrs_dim1];
	hwidth[j] = hwidts[j + *sbrgns * hwidts_dim1];
/* L10: */
    }

/*    If DVFLAG = 1, we will remove the region */
/*    with greatest estimated error from the heap. */

    if (*dvflag == 1) {
	--(*sbrgns);
	subrgn = 1;
L20:
	subtmp = subrgn << 1;
	if (subtmp <= *sbrgns) {
	    if (subtmp != *sbrgns) {

/*   Find max. of left and right child. */

		if (greate[subtmp] < greate[subtmp + 1]) {
		    ++subtmp;
		}
	    }

/*   Compare max.child with parent. */
/*   If parent is max., then done. */

	    if (great < greate[subtmp]) {

/*   Move the values at position subtmp up the heap. */

		greate[subrgn] = greate[subtmp];
		i__1 = *numfun;
		for (j = 1; j <= i__1; ++j) {
		    errors[j + subrgn * errors_dim1] = errors[j + subtmp * 
			    errors_dim1];
		    values[j + subrgn * values_dim1] = values[j + subtmp * 
			    values_dim1];
/* L25: */
		}
		dir[subrgn] = dir[subtmp];
		i__1 = *ndim;
		for (j = 1; j <= i__1; ++j) {
		    centrs[j + subrgn * centrs_dim1] = centrs[j + subtmp * 
			    centrs_dim1];
		    hwidts[j + subrgn * hwidts_dim1] = hwidts[j + subtmp * 
			    hwidts_dim1];
/* L30: */
		}
		subrgn = subtmp;
		goto L20;
	    }
	}
    } else if (*dvflag == 2) {

/*   If DVFLAG = 2, then insert new region in the heap. */

	subrgn = *sbrgns;
L40:
	subtmp = subrgn / 2;
	if (subtmp >= 1) {

/*   Compare max.child with parent. */
/*   If parent is max, then done. */

	    if (great > greate[subtmp]) {

/*   Move the values at position subtmp down the heap. */

		greate[subrgn] = greate[subtmp];
		i__1 = *numfun;
		for (j = 1; j <= i__1; ++j) {
		    errors[j + subrgn * errors_dim1] = errors[j + subtmp * 
			    errors_dim1];
		    values[j + subrgn * values_dim1] = values[j + subtmp * 
			    values_dim1];
/* L45: */
		}
		dir[subrgn] = dir[subtmp];
		i__1 = *ndim;
		for (j = 1; j <= i__1; ++j) {
		    centrs[j + subrgn * centrs_dim1] = centrs[j + subtmp * 
			    centrs_dim1];
		    hwidts[j + subrgn * hwidts_dim1] = hwidts[j + subtmp * 
			    hwidts_dim1];
/* L50: */
		}
		subrgn = subtmp;
		goto L40;
	    }
	}
    }

/*    Insert the saved values in their correct places. */

    if (*sbrgns > 0) {
	greate[subrgn] = great;
	i__1 = *numfun;
	for (j = 1; j <= i__1; ++j) {
	    errors[j + subrgn * errors_dim1] = error[j];
	    values[j + subrgn * values_dim1] = value[j];
/* L55: */
	}
	dir[subrgn] = direct;
	i__1 = *ndim;
	for (j = 1; j <= i__1; ++j) {
	    centrs[j + subrgn * centrs_dim1] = center[j];
	    hwidts[j + subrgn * hwidts_dim1] = hwidth[j];
/* L60: */
	}
    }

/* ***END DTRHRE */

    return 0;
} /* dtrhre_ */

#ifdef __cplusplus
	}
#endif
