#!/usr/bin/python
#
# File name:   mklcmo.py
# Date:        2016/05/01 20:53
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np
import random as rn

if __name__ == "__main__":

    nl = 8 - int(sys.argv[1])

    #nx = 36
    #ny = 36
    nx = 8
    ny = 8
    nz = 18
    #nx = 20
    #ny = 20
    #nz = 3
    ti = []
    xi = []
    yi = []
    zi = []

    for iy in range(ny):
        for ix in range(nx):
            tt = 3
            ti.append(tt)
            xi.append(ix+1)
            yi.append(iy+1)
            zi.append(1)

    iii = 0
    for iz in range(1,nz-1):
        for iy in range(ny):
            for ix in range(nx):
                iii = iii + 1
                if iii % 8 < nl:
                    tt = 1
                else:
                    tt = 2
                ti.append(tt)
                xi.append(ix+1)
                yi.append(iy+1)
                zi.append(iz+1)

    for iy in range(ny):
        for ix in range(nx):
            tt = 0
            ti.append(tt)
            xi.append(ix+1)
            yi.append(iy+1)
            zi.append(nz)

    np = len(ti)

    print np
    print nx, ny, nz
    for i in range(np):
        print ti[i], xi[i], yi[i], zi[i]
        #if ti[i] == 1:
        #    a = 'Ca'
        #elif ti[i] == 2:
        #    a = 'N'
        #elif ti[i] == 3:
        #    a = 'F'
        #print a, xi[i], yi[i], zi[i]

# end of mklcmo.py 
