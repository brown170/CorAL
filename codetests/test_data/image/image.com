image_controls{
   # knot_init_scheme = default # or no_data, sampling_thm, user_defined, optimized 
    knot_init_scheme = optimized 
    optimize_knots = true
    constrain_origin = true
    constrain_rmax_zero = false
    qmax = 70
    rmax = 60 
    numcoeffs = 10
    qscale = 70.
    hbt_only = true
    knot_tolerance = 3.e-4
}
