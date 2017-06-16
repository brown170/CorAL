8/1/2012
David Brown

==================================================================
                    Code Tests for CorAL 
==================================================================

All these tests build with the main SConstruct file in the top-level
directory of the project:

    > scons tests=True

They install to the codetests/ directory, along with their source
and all needed data files.


-------------------- unit tests & the like -------------------- 

Integrator tests:

    dcuhre/dtest1 -- FAILS TO BUILD::
    dcuhre/dtest2 -- FAILS TO BUILD::


Linear algebra tests:

    complexmatrixtest -- linear algebra tests or complex N-dimensional matrices
    realmatrixtest -- linear algebra tests or real N-dimensional matrices
    tnt_inversions -- matrix math and GLSQR inverter tests


Random numbers:

    gslrantest -- test GSL's random number generator
    rantest -- test CorAL's random number generator


Special function tests:

    chtester -- Cartesian harmonics
    gslbesstester -- test of GSL's Bessel functions
    gslbug -- test of a GSL bug in its Coulomb WF package
    ylmchtester -- Spherical <--> Cartesian harmonic conversions
    ylmtester -- a Spherical harmonic (Ylm) tester


Tests of miscellaneous utilities:

    dycast_tester -- test various smart-pointer tools used in CorAL
    testBoost -- Lorentz boosts


Wavefunction test codes:

    coulwave/cgtester -- 
    coulwave/cwtester -- FAILS TO BUILD:: function CWincoming_smallr not declared in wavefunction.h
    iwtest -- 
    phaseshift_tester -- ??             \___ which is which?
    sortTheseOut/phseshift_tester -- ?? /
    sortTheseOut/wftestX -- ?? 
    planewavetest -- 
    wfsample -- 


stester -- simple CorAL regression tests: type "stester -help" to find out details

    Valid tests are: -pmap       : Test read/write/assignment of parameterMap
                     -objio      : Read/write various objects
                     -harm       : FAILS:: Initialize a 3d correlation represented in a Spherical 
                                   Harmonic basis
                     -image      : FAILS:: Simple basis spline imaging test in 1d
                     -legimage   : FAILS:: Simple Legendre polynomial imaging test in 1d
                     -legimage3d : FAILS:: Legendre polynomial imaging test in 3d
                     -conv       : FAILS:: Convolute a Gaussian with a kernel

    This code uses all tests in the imageLegendreTest.h, objectIOTest.h, 
    parameterMapTest.h headers.


Obsolete tests:

    besstester -- test of Bessel functions
    coulwave/ptester -- 
    fixUncert -- 
    ylm_tester -- a Spherical harmonic (Ylm) tester
    
    
-------------------- integration tests -------------------- 


Kernels:

    kernelsample -- 
    sortTheseOut/kerneltestX -- FAILS TO BUILD::
    kzeros -- zeros of the K_l=0(q,r) and K_l=2(q,r) pi0-pi0 kernel


Miscellaneous:

    boltztest -- 
    cgctester -- 
    testConvolve -- 


Source functions:

    sourcesample_blast --
    sourcesample_gauss --
    sourcesample_OSCAR -- FAILS TO BUILD:: wrong function prototype for CSourceCalc_OSCAR::CalcS ??


    


