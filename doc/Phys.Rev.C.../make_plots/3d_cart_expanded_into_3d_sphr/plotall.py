import os,glob

for t in glob.glob("source/*.coral"):
    p = t.replace(".coral","_plot.dat")
    if p.find("hermite")!=-1:       mode = "hermite"
    elif p.find("laguerre")!=-1:    mode = "laguerre"
    elif p.find("legendre")!=-1:    mode = "legendre"
    elif p.find("bspline")!=-1:     mode = "bspline"
    elif p.find("chebyshev")!=-1:   mode = "chebyshev"
    else:                           mode = "histogram"
    os.system("../scplot -"+mode+" "+t)
    if os.path.exists("output_plot.dat"):
        os.system("mv output_plot.dat "+p)
    else: print "Generation of "+p+" failed!"

for t in glob.glob("correlation/*.coral"):
    p = t.replace(".coral","_plot.dat")
    mode = "histogram"
    os.system("../scplot -"+mode+" "+t)
    if os.path.exists("output_plot.dat"):
        os.system("mv output_plot.dat "+p)
    else: print "Generation of "+p+" failed!"
