/* dcuhre.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Table of constant values */

static integer c__15 = 15;
static integer c__1 = 1;

/* Subroutine */ int dcuhre_(integer *ndim, integer *numfun, doublereal *a, 
	doublereal *b, integer *minpts, integer *maxpts, U_fp funsub, 
	doublereal *epsabs, doublereal *epsrel, integer *key, integer *nw, 
	integer *restar, doublereal *result, doublereal *abserr, integer *
	neval, integer *ifail, doublereal *work)
{
    static integer i1, i2, i3, i4, i5, i6, i7, i8, k1, k2, k3, k4, k5, k6, k7,
	     k8, num, keyf, lenw, nsub;
    static doublereal work2[648];
    extern /* Subroutine */ int dadhre_(integer *, integer *, integer *, 
	    doublereal *, doublereal *, integer *, integer *, U_fp, 
	    doublereal *, doublereal *, integer *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *), dchhre_(integer *, integer *, integer *, integer *,
	     doublereal *, doublereal *, integer *, integer *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *);
    static integer minsub, maxsub, wtleng, wrksub;

/* ***BEGIN PROLOGUE DCUHRE */
/* ***DATE WRITTEN   900116   (YYMMDD) */
/* ***REVISION DATE  900116   (YYMMDD) */
/* ***CATEGORY NO. H2B1A1 */
/* ***AUTHOR */
/*            Jarle Berntsen, The Computing Centre, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, Norway */
/*            Phone..  47-5-544055 */
/*            Email..  jarle@eik.ii.uib.no */
/*            Terje O. Espelid, Department of Informatics, */
/*            University of Bergen, Thormohlens gt. 55, */
/*            N-5008 Bergen, Norway */
/*            Phone..  47-5-544180 */
/*            Email..  terje@eik.ii.uib.no */
/*            Alan Genz, Computer Science Department, Washington State */
/*            University, Pullman, WA 99163-1210, USA */
/*            Phone.. 509-335-2131 */
/*            Email..  acg@cs2.cs.wsu.edu */
/* ***KEYWORDS automatic multidimensional integrator, */
/*            n-dimensional hyper-rectangles, */
/*            general purpose, global adaptive */
/* ***PURPOSE  The routine calculates an approximation to a given */
/*            vector of definite integrals */

/*      B(1) B(2)     B(NDIM) */
/*     I    I    ... I       (F ,F ,...,F      ) DX(NDIM)...DX(2)DX(1), */
/*      A(1) A(2)     A(NDIM)  1  2      NUMFUN */

/*       where F = F (X ,X ,...,X    ), I = 1,2,...,NUMFUN. */
/*              I   I  1  2      NDIM */

/*            hopefully satisfying for each component of I the following */
/*            claim for accuracy: */
/*            ABS(I(K)-RESULT(K)).LE.MAX(EPSABS,EPSREL*ABS(I(K))) */
/* ***DESCRIPTION Computation of integrals over hyper-rectangular */
/*            regions. */
/*            DCUHRE is a driver for the integration routine */
/*            DADHRE, which repeatedly subdivides the region */
/*            of integration and estimates the integrals and the */
/*            errors over the subregions with greatest */
/*            estimated errors until the error request */
/*            is met or MAXPTS function evaluations have been used. */

/*            For NDIM = 2 the default integration rule is of */
/*            degree 13 and uses 65 evaluation points. */
/*            For NDIM = 3 the default integration rule is of */
/*            degree 11 and uses 127 evaluation points. */
/*            For NDIM greater then 3 the default integration rule */
/*            is of degree 9 and uses NUM evaluation points where */
/*              NUM = 1 + 4*2*NDIM + 2*NDIM*(NDIM-1) + 4*NDIM*(NDIM-1) + */
/*                    4*NDIM*(NDIM-1)*(NDIM-2)/3 + 2**NDIM */
/*            The degree 9 rule may also be applied for NDIM = 2 */
/*            and NDIM = 3. */
/*            A rule of degree 7 is available in all dimensions. */
/*            The number of evaluation */
/*            points used by the degree 7 rule is */
/*              NUM = 1 + 3*2*NDIM + 2*NDIM*(NDIM-1) + 2**NDIM */

/*            When DCUHRE computes estimates to a vector of */
/*            integrals, all components of the vector are given */
/*            the same treatment. That is, I(F ) and I(F ) for */
/*                                            J         K */
/*            J not equal to K, are estimated with the same */
/*            subdivision of the region of integration. */
/*            For integrals with enough similarity, we may save */
/*            time by applying DCUHRE to all integrands in one call. */
/*            For integrals that vary continuously as functions of */
/*            some parameter, the estimates produced by DCUHRE will */
/*            also vary continuously when the same subdivision is */
/*            applied to all components. This will generally not be */
/*            the case when the different components are given */
/*            separate treatment. */

/*            On the other hand this feature should be used with */
/*            caution when the different components of the integrals */
/*            require clearly different subdivisions. */

/*   ON ENTRY */

/*     NDIM   Integer. */
/*            Number of variables. 1 < NDIM <=  15. */
/*     NUMFUN Integer. */
/*            Number of components of the integral. */
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
/*            For 3 < NDIM < 13 the minimum values for MAXPTS are: */
/*             NDIM =    4   5   6    7    8    9    10   11    12 */
/*            KEY = 3:  459 819 1359 2151 3315 5067 7815 12351 20235 */
/*            KEY = 4:  195 309  483  765 1251 2133 3795  7005 13299 */
/*     FUNSUB Externally declared subroutine for computing */
/*            all components of the integrand at the given */
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
/*            In this case the only parameters for DCUHRE that may */
/*            be changed (with respect to the previous call of DCUHRE) */
/*            are MINPTS, MAXPTS, EPSABS, EPSREL and RESTAR. */

/*   ON RETURN */

/*     RESULT Real array of dimension NUMFUN. */
/*            Approximations to all components of the integral. */
/*     ABSERR Real array of dimension NUMFUN. */
/*            Estimates of absolute errors. */
/*     NEVAL  Integer. */
/*            Number of function evaluations used by DCUHRE. */
/*     IFAIL  Integer. */
/*            IFAIL = 0 for normal exit, when ABSERR(K) <=  EPSABS or */
/*              ABSERR(K) <=  ABS(RESULT(K))*EPSREL with MAXPTS or less */
/*              function evaluations for all values of K, */
/*              1 <= K <= NUMFUN . */
/*            IFAIL = 1 if MAXPTS was too small for DCUHRE */
/*              to obtain the required accuracy. In this case DCUHRE */
/*              returns values of RESULT with estimated absolute */
/*              errors ABSERR. */
/*            IFAIL = 2 if KEY is less than 0 or KEY greater than 4. */
/*            IFAIL = 3 if NDIM is less than 2 or NDIM greater than 15. */
/*            IFAIL = 4 if KEY = 1 and NDIM not equal to 2. */
/*            IFAIL = 5 if KEY = 2 and NDIM not equal to 3. */
/*            IFAIL = 6 if NUMFUN is less than 1. */
/*            IFAIL = 7 if volume of region of integration is zero. */
/*            IFAIL = 8 if MAXPTS is less than 3*NUM. */
/*            IFAIL = 9 if MAXPTS is less than MINPTS. */
/*            IFAIL = 10 if EPSABS < 0 and EPSREL < 0. */
/*            IFAIL = 11 if NW is too small. */
/*            IFAIL = 12 if unlegal RESTAR. */
/*     WORK   Real array of dimension NW. */
/*            Used as working storage. */
/*            WORK(NW) = NSUB, the number of subregions in the data */
/*            structure. */
/*            Let WRKSUB=(NW-1-17*NUMFUN*MDIV)/(2*NDIM+2*NUMFUN+2) */
/*            WORK(1),...,WORK(NUMFUN*WRKSUB) contain */
/*              the estimated components of the integrals over the */
/*              subregions. */
/*            WORK(NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB) contain */
/*              the estimated errors over the subregions. */
/*            WORK(2*NUMFUN*WRKSUB+1),...,WORK(2*NUMFUN*WRKSUB+NDIM* */
/*              WRKSUB) contain the centers of the subregions. */
/*            WORK(2*NUMFUN*WRKSUB+NDIM*WRKSUB+1),...,WORK((2*NUMFUN+ */
/*              NDIM)*WRKSUB+NDIM*WRKSUB) contain subregion half widths. */
/*            WORK(2*NUMFUN*WRKSUB+2*NDIM*WRKSUB+1),...,WORK(2*NUMFUN* */
/*              WRKSUB+2*NDIM*WRKSUB+WRKSUB) contain the greatest errors */
/*              in each subregion. */
/*            WORK((2*NUMFUN+2*NDIM+1)*WRKSUB+1),...,WORK((2*NUMFUN+ */
/*              2*NDIM+1)*WRKSUB+WRKSUB) contain the direction of */
/*              subdivision in each subregion. */
/*            WORK(2*(NDIM+NUMFUN+1)*WRKSUB),...,WORK(2*(NDIM+NUMFUN+1)* */
/*              WRKSUB+ 17*MDIV*NUMFUN) is used as temporary */
/*              storage in DADHRE. */


/*        DCUHRE Example Test Program */


/*   DTEST1 is a simple test driver for DCUHRE. */

/*   Output produced on a SUN 3/50. */

/*       DCUHRE TEST RESULTS */

/*    FTEST CALLS = 3549, IFAIL =  0 */
/*   N   ESTIMATED ERROR   INTEGRAL */
/*   1     0.00000010     0.13850818 */
/*   2     0.00000013     0.06369469 */
/*   3     0.00000874     0.05861748 */
/*   4     0.00000021     0.05407034 */
/*   5     0.00000019     0.05005614 */
/*   6     0.00000009     0.04654608 */

/*     PROGRAM DTEST1 */
/*     EXTERNAL FTEST */
/*     INTEGER KEY, N, NF, NDIM, MINCLS, MAXCLS, IFAIL, NEVAL, NW */
/*     PARAMETER (NDIM = 5, NW = 5000, NF = NDIM+1) */
/*     DOUBLE PRECISION A(NDIM), B(NDIM), WRKSTR(NW) */
/*     DOUBLE PRECISION ABSEST(NF), FINEST(NF), ABSREQ, RELREQ */
/*     DO 10 N = 1,NDIM */
/*        A(N) = 0 */
/*        B(N) = 1 */
/*  10 CONTINUE */
/*     MINCLS = 0 */
/*     MAXCLS = 10000 */
/*     KEY = 0 */
/*     ABSREQ = 0 */
/*     RELREQ = 1E-3 */
/*     CALL DCUHRE(NDIM, NF, A, B, MINCLS, MAXCLS, FTEST, ABSREQ, RELREQ, */
/*    * KEY, NW, 0, FINEST, ABSEST, NEVAL, IFAIL, WRKSTR) */
/*     PRINT 9999, NEVAL, IFAIL */
/* 9999 FORMAT (8X, 'DCUHRE TEST RESULTS', //'     FTEST CALLS = ', I4, */
/*    * ', IFAIL = ', I2, /'    N   ESTIMATED ERROR   INTEGRAL') */
/*     DO 20 N = 1,NF */
/*        PRINT 9998, N, ABSEST(N), FINEST(N) */
/* 9998    FORMAT (3X, I2, 2F15.8) */
/*  20 CONTINUE */
/*     END */
/*     SUBROUTINE FTEST(NDIM, Z, NFUN, F) */
/*     INTEGER N, NDIM, NFUN */
/*     DOUBLE PRECISION Z(NDIM), F(NFUN), SUM */
/*     SUM = 0 */
/*     DO 10 N = 1,NDIM */
/*        SUM = SUM + N*Z(N)**2 */
/*  10 CONTINUE */
/*     F(1) = EXP(-SUM/2) */
/*     DO 20 N = 1,NDIM */
/*        F(N+1) = Z(N)*F(1) */
/*  20 CONTINUE */
/*     END */

/* ***LONG DESCRIPTION */

/*   The information for each subregion is contained in the */
/*   data structure WORK. */
/*   When passed on to DADHRE, WORK is split into eight */
/*   arrays VALUES, ERRORS, CENTRS, HWIDTS, GREATE, DIR, */
/*   OLDRES and WORK. */
/*   VALUES contains the estimated values of the integrals. */
/*   ERRORS contains the estimated errors. */
/*   CENTRS contains the centers of the subregions. */
/*   HWIDTS contains the half widths of the subregions. */
/*   GREATE contains the greatest estimated error for each subregion. */
/*   DIR    contains the directions for further subdivision. */
/*   OLDRES and WORK are used as work arrays in DADHRE. */

/*   The data structures for the subregions are in DADHRE organized */
/*   as a heap, and the size of GREATE(I) defines the position of */
/*   region I in the heap. The heap is maintained by the program */
/*   DTRHRE. */

/*   The subroutine DADHRE is written for efficient execution on shared */
/*   memory parallel computer. On a computer with NPROC processors we wil */
/*   in each subdivision step divide MDIV regions, where MDIV is */
/*   chosen such that MOD(2*MDIV,NPROC) = 0, in totally 2*MDIV new region */
/*   Each processor will then compute estimates of the integrals and erro */
/*   over 2*MDIV/NPROC subregions in each subdivision step. */
/*   The subroutine for estimating the integral and the error over */
/*   each subregion, DRLHRE, uses WORK2 as a work array. */
/*   We must make sure that each processor writes its results to */
/*   separate parts of the memory, and therefore the sizes of WORK and */
/*   WORK2 are functions of MDIV. */
/*   In order to achieve parallel processing of subregions, compiler */
/*   directives should be placed in front of the DO 200 */
/*   loop in DADHRE on machines like Alliant and CRAY. */

/* ***REFERENCES */
/*   J.Berntsen, T.O.Espelid and A.Genz, An Adaptive Algorithm */
/*   for the Approximate Calculation of Multiple Integrals, */
/*   To be published. */

/*   J.Berntsen, T.O.Espelid and A.Genz, DCUHRE: An Adaptive */
/*   Multidimensional Integration Routine for a Vector of */
/*   Integrals, To be published. */

/* ***ROUTINES CALLED DCHHRE,DADHRE */
/* ***END PROLOGUE DCUHRE */

/*   Global variables. */


/*   Local variables. */

/*   MDIV   Integer. */
/*          MDIV is the number of subregions that are divided in */
/*          each subdivision step in DADHRE. */
/*          MDIV is chosen default to 1. */
/*          For efficient execution on parallel computers */
/*          with NPROC processors MDIV should be set equal to */
/*          the smallest integer such that MOD(2*MDIV,NPROC) = 0. */
/*   MAXDIM Integer. */
/*          The maximum allowed value of NDIM. */
/*   MAXWT  Integer. The maximum number of weights used by the */
/*          integration rule. */
/*   WTLENG Integer. */
/*          The number of generators used by the selected rule. */
/*   WORK2  Real work space. The length */
/*          depends on the parameters MDIV,MAXDIM and MAXWT. */
/*   MAXSUB Integer. */
/*          The maximum allowed number of subdivisions */
/*          for the given values of KEY, NDIM and MAXPTS. */
/*   MINSUB Integer. */
/*          The minimum allowed number of subregions for the given */
/*          values of MINPTS, KEY and NDIM. */
/*   WRKSUB Integer. */
/*          The maximum allowed number of subregions as a function */
/*          of NW, NUMFUN, NDIM and MDIV. This determines the length */
/*          of the main work arrays. */
/*   NUM    Integer. The number of integrand evaluations needed */
/*          over each subregion. */

/* cccccccccccccccccccccccccc */

/* dave's testing lines */

/*      write(*,*) 'In DCUHRE:  ndim=',NDIM,' numfun=',numfun */
/*      write(*,*) 'In DCUHRE:  minpts=',MINPTS,' maxpts=',MAXPTS */
/*      write(*,*) 'In DCUHRE:  A=',(A(i),i=1,3) */
/*      write(*,*) 'In DCUHRE:  B=',(B(i),i=1,3) */
/*      write(*,*) 'In DCUHRE:  EPSABS=',EPSABS,' EPSREL=',EPSREL */
/*      write(*,*) 'In DCUHRE:  key=',key,' nw=',nw */
/*      call FUNSUB(ndim,b,numfun,result) */
/*      write(*,*) 'In DCUHRE:  FUNSUB= ',(result(i),i=1,numfun) */
/* cccccccccccccccccccccccccc */

/* ***FIRST EXECUTABLE STATEMENT DCUHRE */

/*   Compute NUM, WTLENG, MAXSUB and MINSUB, */
/*   and check the input parameters. */

    /* Parameter adjustments */
    --b;
    --a;
    --abserr;
    --result;
    --work;

    /* Function Body */
    dchhre_(&c__15, ndim, numfun, &c__1, &a[1], &b[1], minpts, maxpts, epsabs,
	     epsrel, key, nw, restar, &num, &maxsub, &minsub, &keyf, ifail, &
	    wtleng);
    wrksub = (*nw - 1 - *numfun * 17) / ((*ndim << 1) + (*numfun << 1) + 2);
    if (*ifail != 0) {
	goto L999;
    }

/*   Split up the work space. */

    i1 = 1;
    i2 = i1 + wrksub * *numfun;
    i3 = i2 + wrksub * *numfun;
    i4 = i3 + wrksub * *ndim;
    i5 = i4 + wrksub * *ndim;
    i6 = i5 + wrksub;
    i7 = i6 + wrksub;
    i8 = i7 + *numfun;
    k1 = 1;
    k2 = k1 + (wtleng << 1) * *ndim;
    k3 = k2 + wtleng * 5;
    k4 = k3 + wtleng;
    k5 = k4 + *ndim;
    k6 = k5 + *ndim;
    k7 = k6 + (*ndim << 1);
    k8 = k7 + wtleng * 3;

/*   On restart runs the number of subregions from the */
/*   previous call is assigned to NSUB. */

    if (*restar == 1) {
	nsub = (integer) work[*nw];
    }

/*   Compute the size of the temporary work space needed in DADHRE. */

    lenw = *numfun << 4;
    dadhre_(ndim, numfun, &c__1, &a[1], &b[1], &minsub, &maxsub, (U_fp)funsub,
	     epsabs, epsrel, &keyf, restar, &num, &lenw, &wtleng, &result[1], 
	    &abserr[1], neval, &nsub, ifail, &work[i1], &work[i2], &work[i3], 
	    &work[i4], &work[i5], &work[i6], &work[i7], &work[i8], &work2[k1 
	    - 1], &work2[k2 - 1], &work2[k3 - 1], &work2[k4 - 1], &work2[k5 - 
	    1], &work2[k6 - 1], &work2[k7 - 1], &work2[k8 - 1]);
    work[*nw] = (doublereal) nsub;
L999:
    return 0;

/* ***END DCUHRE */

} /* dcuhre_ */

#ifdef __cplusplus
	}
#endif
