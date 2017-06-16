/* d132re.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Subroutine */ int d132re_(integer *wtleng, doublereal *w, doublereal *g, 
	doublereal *errcof, doublereal *rulpts)
{
    /* Initialized data */

    static doublereal dim2g[16] = { .2517129343453109,.7013933644534266,
	    .9590960631619962,.9956010478552127,.5,.1594544658297559,
	    .3808991135940188,.6582769255267192,.8761473165029315,
	    .998243184053198,.9790222658168462,.6492284325645389,
	    .8727421201131239,.3582614645881228,.5666666666666666,
	    .2077777777777778 };
    static doublereal dim2w[70]	/* was [14][5] */ = { .0337969236013446,
	    .09508589607597761,.1176006468056962,.0265777458632695,
	    .0170144177020064,0.,.0162659309863741,.1344892658526199,
	    .1328032165460149,.0563747476999187,.0039082790813105,
	    .0301279877743215,.1030873234689166,.0625,.3213775489050763,
	    -.1767341636743844,.07347600537466072,-.03638022004364754,
	    .02125297922098712,.1460984204026913,.01747613286152099,
	    .1444954045641582,1.307687976001325e-4,5.380992313941161e-4,
	    1.042259576889814e-4,-.001401152865045733,.008041788181514763,
	    -.1420416552759383,.3372900883288987,-.1644903060344491,
	    .07707849911634622,-.0380447835850631,.02223559940380806,
	    .1480693879765931,4.467143702185814e-6,.150894476707413,
	    3.647200107516215e-5,5.77719899901388e-4,1.041757313688177e-4,
	    -.001452822267047819,.008338339968783705,-.147279632923196,
	    -.8264123822525677,.306583861409436,.002389292538329435,
	    -.1343024157997222,.088333668405339,0.,9.786283074168292e-4,
	    -.1319227889147519,.00799001220015063,.003391747079760626,
	    .002294915718283264,-.01358584986119197,.04025866859057809,
	    .003760268580063992,.6539094339575232,-.2041614154424632,
	    -.174698151579499,.03937939671417803,.006974520545933992,0.,
	    .006667702171778258,.05512960621544304,.05443846381278607,
	    .02310903863953934,.01506937747477189,-.0605702164890189,
	    .04225737654686337,.02561989142123099 };

    static integer i__, j;

/* ***BEGIN PROLOGUE D132RE */
/* ***AUTHOR   Jarle Berntsen, EDB-senteret, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, NORWAY */
/* ***PURPOSE D132RE computes abscissas and weights of a 2 dimensional */
/*            integration rule of degree 13. */
/*            Two null rules of degree 11, one null rule of degree 9 */
/*            and one null rule of degree 7 to be used in error */
/*            estimation are also computed. */
/* ***DESCRIPTION D132RE will select the correct values of the abscissas */
/*            and corresponding weights for different */
/*            integration rules and null rules and assign them to */
/*            G and W. The heuristic error coefficients ERRCOF */
/*            will also be assigned. */


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
/*            The number of points produced by each generator. */
/* ***REFERENCES S.Eriksen, */
/*              Thesis of the degree cand.scient, Dept. of Informatics, */
/*              Univ. of Bergen,Norway, 1984. */

/* ***ROUTINES CALLED-NONE */
/* ***END PROLOGUE D132RE */

/*   Global variables */


/*   Local variables. */


    /* Parameter adjustments */
    --rulpts;
    g -= 3;
    w -= 6;
    --errcof;

    /* Function Body */






/* ***FIRST EXECUTABLE STATEMENT D132RE */

/*   Assign values to W. */

    for (i__ = 1; i__ <= 14; ++i__) {
	for (j = 1; j <= 5; ++j) {
	    w[j + i__ * 5] = dim2w[i__ + j * 14 - 15];
/* L10: */
	}
    }

/*   Assign values to G. */

    for (i__ = 1; i__ <= 2; ++i__) {
	for (j = 1; j <= 14; ++j) {
	    g[i__ + (j << 1)] = 0.;
/* L20: */
	}
    }
    g[5] = dim2g[0];
    g[7] = dim2g[1];
    g[9] = dim2g[2];
    g[11] = dim2g[3];
    g[13] = dim2g[4];
    g[15] = dim2g[5];
    g[16] = g[15];
    g[17] = dim2g[6];
    g[18] = g[17];
    g[19] = dim2g[7];
    g[20] = g[19];
    g[21] = dim2g[8];
    g[22] = g[21];
    g[23] = dim2g[9];
    g[24] = g[23];
    g[25] = dim2g[10];
    g[26] = dim2g[11];
    g[27] = dim2g[12];
    g[28] = dim2g[13];
    g[29] = dim2g[14];
    g[30] = dim2g[15];

/*   Assign values to RULPTS. */

    rulpts[1] = 1.;
    for (i__ = 2; i__ <= 11; ++i__) {
	rulpts[i__] = 4.;
/* L30: */
    }
    rulpts[12] = 8.;
    rulpts[13] = 8.;
    rulpts[14] = 8.;

/*   Assign values to ERRCOF. */

    errcof[1] = 10.;
    errcof[2] = 10.;
    errcof[3] = (float)1.;
    errcof[4] = (float)5.;
    errcof[5] = (float).5;
    errcof[6] = (float).25;

/* ***END D132RE */

    return 0;
} /* d132re_ */

#ifdef __cplusplus
	}
#endif
