#!/usr/bin/python
#
# File name:   sim_lg.py
# Date:        2016/04/19 22:55
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

    rn.seed(145146167)

# model
    nt = 5
    # valence
    vv = np.zeros((nt), dtype=float)
    vv[0] = -2.0 # O(-II)
    vv[1] =  3.0 # Mn(III)
    vv[2] =  4.0 # Mn(IV)
    vv[3] =  2.0 # Ca(II)
    vv[4] =  3.0 # La(III)

    # interactions
    ee = np.zeros((nt,nt), dtype=float)
    ee[0,1] = -1.0 # O - Mn(III)
    ee[1,0] = ee[0,1]
    ee[0,2] = -1.0 # O - Mn(IV)
    ee[2,0] = ee[0,2]
    ee[0,3] = -1.0 # O - Ca(II)
    ee[3,0] = ee[0,3]
    ee[0,4] = -1.2 # O - La(III)
    ee[4,0] = ee[0,4]

    ee[1,2] = -1.0 # Mn(III) - Mn(IV) resonance

    ee[0,0] = 1.0 # O-O Coulombic repulsion
    ee[3,3] = 1.0 # Ca-Ca Coulombic repulsion
    ee[4,4] = 1.5 # La-La Coulombic repulsion
    ee[3,4] = (ee[3,3]*ee[4,4])**0.5 # Ca-La Coulombic repulsion
    ee[4,3] = ee[3,4]

    # interaction - configurations
    # Mn-O (1s), Mn-Mn (2s)
    # Ca-O (1e), Ca-Ca (2s)
    # Om-Mn (1sp), Om-Ca (1ed)
    # Oc-Ca (1ep), Om-Ca (1e)
    
# configuration
    # read initial configuration
    f = open(sys.argv[1], 'r')

    np = int(re.findall('\S+', f.readline())[0])
    xi = np.zeros((np), dtype=int)
    yi = np.zeros((np), dtype=int)
    zi = np.zeros((np), dtype=int)

    nx, ny, nz = map(int, re.findall('\S+', f.readline())[0:3])
    box = np.zeros((nx,ny,nz), dtype=int)

    for i in range(np):
        ti[i], xi[i], yi[i], zi[i] = map(int, re.findall('\S+', f.readline())[0:4])

    f.close()

    # initial energy
    for i in range(np):
        ti = t[i]
        #neighbors
        #Mn - sides(6)
        #Ca,La - edges(12)
        #O-sides(2 out of 4 Mn),edges(4 in plane out of 12)
        
        for j in range(np):
    
    exit()

    # nvt (uvt) simulation cycle 
    for it in range(niter):
        i = rn.randint(npm) # select a movable particle

        move(i)
        xi = x[i]
        yi = y[i]
        zi = z[i]
        ti = t[i]
        uene(i)

        # trial move (depends on the type)
        # Mn - switch oxidation state
        # O - detach (GCMC), move on surface
        # Ca - detach (GCMC), move on surface, switch with La
        # La - detach (GCMC), move on surface, switch with Ca


# end of sim_lg.py 