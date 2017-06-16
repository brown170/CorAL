/* d113re.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int d113re_(integer *wtleng, doublereal *w, doublereal *g, 
	doublereal *errcof, doublereal *rulpts)
{
    /* Initialized data */

    static doublereal dim3g[14] = { .19,.5,.75,.8,.9949999999999999,
	    .99873449983514,.7793703685672423,.9999698993088767,
	    .7902637224771788,.4403396687650737,.4378478459006862,
	    .9549373822794593,.9661093133630748,.4577105877763134 };
    static doublereal dim3w[65]	/* was [13][5] */ = { .007923078151105734,
	    .0679717739278808,.001086986538805825,.1838633662212829,
	    .03362119777829031,.01013751123334062,.001687648683985235,
	    .1346468564512807,.001750145884600386,.07752336383837454,
	    .2461864902770251,.06797944868483039,.01419962823300713,
	    1.715006248224684,-.3755893815889209,.1488632145140549,
	    -.2497046640620823,.1792501419135204,.00344612675897389,
	    -.005140483185555825,.006536017839876425,-6.5134549392297e-4,
	    -.006304672433547204,.01266959399788263,-.005454241018647931,
	    .004826995274768427,1.936014978949526,-.3673449403754268,
	    .02929778657898176,-.1151883520260315,.05086658220872218,
	    .04453911087786469,-.022878282571259,.02908926216345833,
	    -.002898884350669207,-.02805963413307495,.05638741361145884,
	    -.02427469611942451,.02148307034182882,.517082819560576,
	    .01445269144914044,-.3601489663995932,.3628307003418485,
	    .007148802650872729,-.09222852896022966,.01719339732471725,
	    -.102141653746035,-.007504397861080493,.01648362537726711,
	    .05234610158469334,.01445432331613066,.003019236275367777,
	    2.05440450381852,.0137775998849012,-.576806291790441,
	    .03726835047700328,.006814878939777219,.05723169733851849,
	    -.04493018743811285,.02729236573866348,3.54747395055699e-4,
	    .01571366799739551,.04990099219278567,.0137791555266677,
	    .002878206423099872 };

    static integer i__, j;

/* ***BEGIN PROLOGUE D113RE */
/* ***AUTHOR   Jarle Berntsen, EDB-senteret, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, NORWAY */
/* ***PURPOSE D113RE computes abscissas and weights of a 3 dimensional */
/*            integration rule of degree 11. */
/*            Two null rules of degree 9, one null rule of degree 7 */
/*            and one null rule of degree 5 to be used in error */
/*            estimation are also computed. */
/* ***DESCRIPTION D113RE will select the correct values of the abscissas */
/*            and corresponding weights for different */
/*            integration rules and null rules and assign them to G */
/*            and W. */
/*            The heuristic error coefficients ERRCOF */
/*            will also be computed. */


/*   ON ENTRY */

/*     WTLENG Integer. */
/*            The number of weights in each of the rules. */

/*   ON RETURN */

/*     W      Real array of dimension (5,WTLENG). */
/*            The weights for the basic and null rules. */
/*            W(1,1),...,W(1,WTLENG) are weights for the basic rule. */
/*            W(I,1),...,W(I,WTLENG), for I > 1 are null rule weights. */
/*     G      Real array of dimension (NDIM,WTLENG). */
/*            The fully symmetric sum generators for the rules. */
/*            G(1,J),...,G(NDIM,J) are the generators for the points */
/*            associated with the the Jth weights. */
/*     ERRCOF Real array of dimension 6. */
/*            Heuristic error coefficients that are used in the */
/*            error estimation in BASRUL. */
/*     RULPTS Real array of dimension WTLENG. */
/*            The number of points used by each generator. */

/* ***REFERENCES  J.Berntsen, Cautious adaptive numerical integration */
/*               over the 3-cube, Reports in Informatics 17, Dept. of */
/*               Inf.,Univ. of Bergen, Norway, 1985. */
/*               J.Berntsen and T.O.Espelid, On the construction of */
/*               higher degree three-dimensional embedded integration */
/*               rules, SIAM J. Numer. Anal.,Vol. 25,No. 1, pp.222-234, */
/*               1988. */
/* ***ROUTINES CALLED-NONE */
/* ***END PROLOGUE D113RE */

/*   Global variables. */


/*   Local variables. */


    /* Parameter adjustments */
    --rulpts;
    g -= 4;
    w -= 6;
    --errcof;

    /* Function Body */






/* ***FIRST EXECUTABLE STATEMENT D113RE */

/*   Assign values to W. */

    for (i__ = 1; i__ <= 13; ++i__) {
	for (j = 1; j <= 5; ++j) {
	    w[j + i__ * 5] = dim3w[i__ + j * 13 - 14];
/* L10: */
	}
    }

/*   Assign values to G. */

    for (i__ = 1; i__ <= 3; ++i__) {
	for (j = 1; j <= 13; ++j) {
	    g[i__ + j * 3] = 0.;
/* L20: */
	}
    }
    g[7] = dim3g[0];
    g[10] = dim3g[1];
    g[13] = dim3g[2];
    g[16] = dim3g[3];
    g[19] = dim3g[4];
    g[22] = dim3g[5];
    g[23] = g[22];
    g[25] = dim3g[6];
    g[26] = g[25];
    g[28] = dim3g[7];
    g[29] = g[28];
    g[30] = g[28];
    g[31] = dim3g[8];
    g[32] = g[31];
    g[33] = g[31];
    g[34] = dim3g[9];
    g[35] = g[34];
    g[36] = g[34];
    g[37] = dim3g[11];
    g[38] = dim3g[10];
    g[39] = g[38];
    g[40] = dim3g[12];
    g[41] = g[40];
    g[42] = dim3g[13];

/*   Assign values to RULPTS. */

    rulpts[1] = 1.;
    rulpts[2] = 6.;
    rulpts[3] = 6.;
    rulpts[4] = 6.;
    rulpts[5] = 6.;
    rulpts[6] = 6.;
    rulpts[7] = 12.;
    rulpts[8] = 12.;
    rulpts[9] = 8.;
    rulpts[10] = 8.;
    rulpts[11] = 8.;
    rulpts[12] = 24.;
    rulpts[13] = 24.;

/*   Assign values to ERRCOF. */

    errcof[1] = 4.;
    errcof[2] = (float)4.;
    errcof[3] = (float).5;
    errcof[4] = (float)3.;
    errcof[5] = (float).5;
    errcof[6] = (float).25;

/* ***END D113RE */

    return 0;
} /* d113re_ */

#ifdef __cplusplus
	}
#endif
