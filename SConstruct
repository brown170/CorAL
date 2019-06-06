# <<BEGIN-copyright>>
#
#                 The GNU General Public License (GPL) Version 2, June 1991
#
# Copyright (c) 2013, Lawrence Livermore National Security, LLC. Produced at the Lawrence
# Livermore National Laboratory. Written by Ron Soltz (soltz1@llnl.gov), David A. Brown
# (dbrown@bnl.gov) and Scott Pratt (pratts@pa.msu.edu).
#
# CODE-CODE-643336 All rights reserved.
#
# This file is part of CorAL, Version: 1.17.
#
# Please see the file LICENSE.TXT in the main directory of this source code distribution.
#
# This program is free software; you can redistribute it and/or modify it under the terms of
# the GNU General Public License (as published by the Free Software Foundation) version 2,
# dated June 1991.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the terms and conditions of the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program;
# if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA 02111-1307 USA
#
# <<END-copyright>>
import sys,os

if "CORAL_GSLPATH" in os.environ: gslpath = os.environ['CORAL_GSLPATH']
else:                             gslpath = "/opt/local/"
if "CORAL_X11PATH" in os.environ: x11path = os.environ['CORAL_X11PATH']
else:                             x11path = "/opt/local/"
if "CORAL_HOME" in os.environ: coralhome = os.environ['CORAL_HOME']
else:                          coralhome = os.getcwd()
if "CORAL_CCFLAGS" in os.environ: ccflags = os.environ['CORAL_CCFLAGS']
else:                             ccflags = ""
ccflags += " -Df2cFortran -ansi"

opts = Variables()
opts.Add( BoolVariable( 'debug',    'Set to true to build for debugging', 0 ) )
opts.Add( BoolVariable( 'X11',      'Build libraries dependent on X11', 1 ) )
opts.Add( BoolVariable( 'tests',    'Build regression tests', 0 ) )
opts.Add( BoolVariable( 'doc',      'Build documentation (user guide, api docs)', 0 ) )


# --------------- Default Build Environment ---------------
coralenv = Environment( \
    options = opts,\
    CXX = 'clang',\
    CPPPATH = [ '#include', '#src', gslpath+'include' ],\
    LIBPATH = [ '#lib', gslpath+'lib' ],\
    LIBS    = [ 'gsl', 'gslcblas' ],\
    CCFLAGS = ccflags, \
    LDFLAGS = '',\
    #ENV = \{'PATH':os.environ['PATH']\}\
)

if coralenv['debug']: coralenv['CCFLAGS'] += " -Wall -g -DTNT_BOUNDS_CHECK"

SConscript('src/coralutils/SConscript', exports='coralenv', variant_dir='#build/coralutils', duplicate=0)
SConscript('src/coral/SConscript', exports='coralenv', variant_dir='#build/coral', duplicate=0)

#if coralenv['X11']:
#    coralenv[ 'CPPPATH' ] += [ x11path+'include' ]
#    coralenv[ 'LIBPATH' ] += [ x11path+'lib' ]
#    coralenv[ 'LIBS' ]    += [ 'X11']  # for some reason, adding to this list on MacOSX messes up CheckLib calls later
#
#    # --------------- Build XGraph ---------------
#    SConscript('src/xgraph/SConscript', exports='coralenv', variant_dir='#build/xgraph', duplicate=0)

Help(opts.GenerateHelpText(coralenv))

conf = Configure(coralenv)
#for lib in coralenv['LIBS']:
#    if not conf.CheckLib(lib):
#        print 'Did not find '+lib+' in '+str(coralenv['LIBPATH'])+', exiting!'
#        #Exit(1)
#        coralenv = conf.Finish()

# --------------- Build CorAL & Friends ---------------
SConscript('src/coralutils/SConscript', exports='coralenv', variant_dir='#build/coralutils', duplicate=0)
SConscript('src/coral/SConscript', exports='coralenv', variant_dir='#build/coral', duplicate=0)
#SConscript('src/xgraph/SConscript', exports='coralenv', variant_dir='#build/xgraph', duplicate=0)
SConscript('src/crab/SConscript', exports='coralenv', variant_dir='#build/crab', duplicate=0)
SConscript('src/chum/SConscript', exports='coralenv', variant_dir='#build/chum', duplicate=0)
SConscript('src/diver/SConscript', exports='coralenv', variant_dir='#build/diver', duplicate=0)
SConscript('src/shark/SConscript', exports='coralenv', variant_dir='#build/shark', duplicate=0)
SConscript('src/bfplot/SConscript', exports='coralenv', variant_dir='#build/bfplot', duplicate=0)
SConscript('src/scplot/SConscript', exports='coralenv', variant_dir='#build/scplot', duplicate=0)
SConscript('src/converts2c/SConscript', exports='coralenv', variant_dir='#build/converts2c', duplicate=0)

# --------------- Build regression tests ---------------
if coralenv['tests']:
    SConscript( 'codetests/SConscript', exports = 'coralenv', variant_dir = '#build/tests', duplicate = 0 )

# --------------- Build documentation ---------------
if coralenv['doc']:
    SConscript( 'doc/SConscript', exports = 'coralenv', variant_dir = '#build/doc', duplicate = 0 )
