# CorAL
#--------------------------
# lm = 00 term
#--------------------------
dummyline fred echo {
}
#doesn't work: dummyline fred echo {}
#----------------------------------
# read in the correlation terms 
# we need to process by hand
#----------------------------------
read corr00 stage1_corr_terms_0.dat
#----------------------------------
# image the correlation terms by hand
#----------------------------------
image corr00 sou00{
    rmax = 40.
    bspline_degree=3
    qscale = 40.
    qmax = 60
    constrain_origin
    constrain_rmax_zero
    constrain_rmax_zero_slope
    knot_scheme=sampling_thm 
    dummy {
        4.0
        string
        xx 43
    }
}
#----------------------------------
# write out the source terms 
#----------------------------------
write sou00 stage3_sou_terms_0.dat
#----------------------------------
# unimage the source
#----------------------------------
unimage sou00 corr3_00{
    #override_kernel=nofsi
    bigq = false
    qmax = 60.
    dq = 2.5
    bin_scheme = userdefined
}
#----------------------------------
# write out the unimaged correlation 
#----------------------------------
write corr3_00 stage3_corr_terms_0.dat
#----------------------------------
# generate plots
#----------------------------------
slicerad corr00  plot_corrin_slicerad_00.dat
slicerad corr3_00  plot_corrout_slicerad_00.dat
slicerad sou00 plot_sou_slicerad_00.dat

