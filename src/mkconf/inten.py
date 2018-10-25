#!/usr/bin/python
#
# File name:   inten.py
# Date:        2016/05/24 23:36
# Author:      Lukas Vlcek
#
# Description: 
#

import sys
import string as s
import re
import numpy as np

if __name__ == "__main__":

    dat = []
    f = open(sys.argv[1], 'r')
    for line in f:
        dat.append(map(int, re.findall('\S+', line)))
    f.close()

    na = len(dat)
    dat = np.array(dat)

    ima = np.zeros((108, 104), dtype=float)

    for i in range(na):
        if dat[i,0] == 3:
            nc = 0
            nl = 0
            for j in range(na):
                if dat[j,0] < 3:
                    dx = dat[j,1]-dat[i,1]
                    dy = dat[j,2]-dat[i,2]
                    if abs(dx) == 1 and abs(dy) == 1:
                        if dat[j,0] == 1:
                            nc = nc + 1
                        else:
                            nl = nl + 1

            #print i, nc, nl, nc + nl
            intens = 0

            if nl > 0:
                intens = 1.0
                if nl > 1:
                    intens = intens + 0.3
                if nc > 0:
                    intens = intens + 0.3
                if nl + nc > 3:
                    intens = intens + 0.1
            elif nc > 0:
                intens = 0.6
                if nc > 1:
                    intens = intens + 0.4
                    if nc > 2:
                        intens = intens + 0.1

            k = ( (dat[i,1]/2) + (dat[i,2]/2) ) % 2
            if k == 1:
                intens = -intens

            print dat[i,1]/2, dat[i,2]/2, intens*1.8


# end of inten.py 
