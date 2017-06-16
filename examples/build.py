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
#!/usr/bin/python
import os,sys

global COMPILE_OPTS
global LINK_OPTS
COMPILE_OPTS = []
LINK_OPTS    = []

def buildLib():
    os.chdir("..")
    os.system("scons")
    os.chdir("examples")

def buildTests():
    testCodes = ['stester']
    for testCode in testCodes:
        commandList = [\
            "g++ -g -c -Df2cFortran "+' '.join(COMPILE_OPTS)+" "+testCode+".cc -I../include -o "+testCode+".o",\
            "g++ -g "+testCode+".o "+' '.join(LINK_OPTS)+" -L../lib/ -lm -lstdc++ -lgsl -lgslcblas -lcoral -o "+testCode\
        ]
        for command in commandList: 
            print command
            os.system( command )

if __name__=="__main__":
    global COMPILE_OPTS
    global LINK_OPTS
    if os.environ.has_key('SYS_TYPE') and 'chaos' in os.environ['SYS_TYPE']:
        COMPILE_OPTS.append('-I/usr/gapps/phenix/gsl/include')
        LINK_OPTS.append('-L/usr/gapps/phenix/gsl/lib')
    elif os.uname()[0]=='Darwin':
        COMPILE_OPTS.append('-Wno-long-double')
        COMPILE_OPTS.append('-I/sw/include')
    else: 
        pass
    if not "-nolibs"  in sys.argv: buildLib()
    if not "-notests" in sys.argv: buildTests()
