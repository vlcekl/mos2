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

#0.021551724137931036, 0.5965665236051502, 0.8433179723502304, 0.9575371549893843

if __name__ == "__main__":

    # fraction of Re in Mo
    frac = float(sys.argv[1])

    nx = 64
    ny = 64
    nz = 1
    ti = []
    xi = []
    yi = []
    zi = []

    ntot = 64*64/2
    nre = int(frac*ntot)
    nmo = ntot - nre

    # Set all atom types to Mo
    ntype = np.ones((ntot), dtype=int)

    # Change first nre atoms to Re
    ntype[0:nre] = 2

    # Shuffle the array for random distribution of Re and Mo
    np.random.shuffle(ntype)

    iii = 0
    for iz in range(0,nz):
        for iy in range(ny):
            for ix in range(nx):
                # site on a checker board site
                if (ix + iy) % 2 == 0:
                    tt = ntype[iii]
                    ti.append(tt)
                    xi.append(ix+1)
                    yi.append(iy+1)
                    zi.append(iz+1)
                    iii = iii + 1

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
