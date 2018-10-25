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

    nx = 50
    ny = 50
    nz = 40
    ti = []
    xi = []
    yi = []
    zi = []
    for iz in range(nz):
        for iy in range(ny):
            for ix in range(nx):
                tt = -1
                if iz%2 == 0: #MnO2
                    if iy%4 == 0:
                        if ix%4 == 0:
                            tt = 1 #Mn(III)
                        elif ix%4 == 2:
                            tt = 2 #Mn(IV)
                        else:
                            tt = 0 #O
                    elif iy%4 == 2:
                        if ix%4 == 0:
                            tt = 2 #Mn(IV)
                        elif ix%4 == 2:
                            tt = 1 #Mn(III)
                        else:
                            tt = 0 #O
                    else:
                        if ix%2 == 0:
                            tt = 0 #O
                else: #LaO/CaO
                    if iy%4 == 0:
                        if ix%4 == 0 and iz!=nz-1:
                            tt = 0 
                        elif ix%4 == 2:
                            tt = 0 #O
                    elif iy%4 == 2:
                        if ix%4 == 0:
                            tt = 0 #O on Mn(IV)
                        elif ix%4 == 2 and iz!=nz-1:
                            tt = 0 #O on Mn(III)
                    elif iy%4 == 1 or iy%4 == 3:
                        if ix%2 == 1:
                            if iz!=nz-1:
                                t = rn.randint(0,1)
                            else:
                                t = rn.randint(0,3)

                        if t == 0:
                            tt = 3
                        elif t == 1:
                            tt = 4

                if tt >= 0:
                    ti.append(tt)
                    xi.append(ix)
                    yi.append(iy)
                    zi.append(iz)

    np = len(ti)

    print np
    print nx, ny, nz

    for i in range(np):
        #print ti[i], xi[i], yi[i], zi[i]
        if ti[i] == 0:
            a = 'O'
        elif ti[i] == 1 or ti[i] == 2:
            a = 'Mn'
        elif ti[i] == 3:
            a = 'Ca'
        elif ti[i] == 4:
            a = 'N'
        print a, xi[i], yi[i], zi[i]

# end of mklcmo.py 
