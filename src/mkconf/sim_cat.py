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

def upart_ini(i, box, ee, tti, xt, yt, zt,nx,ny,nz):
    dj = []
    du = []

    # oxygens - attraction
    j = box[(xt+1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt+1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)
 
    j = box[(xt-1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    # cations - repulsion
    j = box[(xt+2)%nx,yt,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-2)%nx,yt,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt+2)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt-2)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]]
        dj.append(j)
        du.append(di)

    return dj, du

def upart_t(i, box, ee, tti, xt, yt, zt,nx,ny,nz):
    dj = []
    du = []

    # old attraction
    k = box[(xi[i]+1)%nx,(yi[i]+1)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    k = box[(xi[i]+1)%nx,(yi[i]-1)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)
 
    k = box[(xi[i]-1)%nx,(yi[i]+1)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    k = box[(xi[i]-1)%nx,(yi[i]-1)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)
    # old repulsion
    k = box[(xi[i]+2)%nx,yi[i],zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    k = box[(xi[i]-2)%nx,yi[i],zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    k = box[xi[i],(yi[i]+2)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    k = box[xi[i],(yi[i]-2)%ny,zi[i]]
    if k >= 0:
        di = -ee[ti[i],ti[k]]
        dj.append(k)
        du.append(di)

    # new attraction
    j = box[(xt+1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt+1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    #print xi[i], yi[i], xt, yt
    # new repulsion
    j = box[(xt+2)%nx,yt,zt]
    if j != i and j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-2)%nx,yt,zt]
    if j != i and j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt+2)%ny,zt]
    if j != i and j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt-2)%ny,zt]
    if j != i and j >= 0:
        di = ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    return dj, du

def upart(i, box, ee, tti, xt, yt, zt,nx,ny,nz):
    dj = []
    du = []

    # oxygens - attraction
    j = box[(xt+1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt+1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)
 
    j = box[(xt-1)%nx,(yt+1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-1)%nx,(yt-1)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    # cations - repulsion
    j = box[(xt+2)%nx,yt,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[(xt-2)%nx,yt,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt+2)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    j = box[xt,(yt-2)%ny,zt]
    if j >= 0:
        di = ee[tti,ti[j]] - ee[ti[i],ti[j]]
        dj.append(j)
        du.append(di)

    return dj, du

def accept(i, j, mi, dj, du, ui):
    # change particle location and type
    if mi == 0:
        box[xi[i], yi[i], zi[i]] = -1
        box[xt, yt, zt] = i
        xi[i] = xt
        yi[i] = yt
    elif mi == 1:
        ti[i] = tti
    elif mi == 2:
        ti[i] = tti
        ti[j] = ttj

    # update particle energies
    ui[i] = ui[i] + np.sum(du)
    #print mi, len(dj), np.sum(du)
    #print du
    for k in range(len(dj)):
        ui[dj[k]] = ui[dj[k]] + du[k]

    return

def cfgout(ti,xi,yi,zi,na, it):
    f = open(str(it)+'.xyz', 'w')
    for i in range(na):
        f.write("{0:d} {1:d} {2:d} {3:d}\n".format(ti[i], xi[i], yi[i], zi[i]))
    f.close()
    return

if __name__ == "__main__":

    rn.seed(145146167)

# model
    nt = 7
    # valence
    vv = np.zeros((nt), dtype=float)
    vv[0] =  0.0 # nothing
    vv[1] =  2.0 # Ca(II)
    vv[2] =  3.0 # La(III)
    vv[3] = -2.0 # O(-II)
    vv[4] =  3.0 # Mn(III)
    vv[5] =  4.0 # Mn(IV)
    vv[6] = -2.0 # O(-II)

    # interactions
    ee = np.zeros((nt+1,nt+1), dtype=float)
    ee[1,3] = -5.0 # O - Ca(II)
    ee[3,1] = ee[1,3]
    ee[2,3] = -6.0 # O - La(III)
    ee[3,2] = ee[2,3]

    ee[1,1] = 4.0 # Ca-Ca Coulombic repulsion
    ee[2,2] = 6.0 # La-La Coulombic repulsion
    ee[1,2] = (ee[1,1]*ee[2,2])**0.5 # Ca-La Coulombic repulsion
    ee[2,1] = ee[1,2]

# Moves
    p_tran = 0.5
    p_mute = 0.0 + p_tran
    p_swap = 0.5 + p_mute

    beta = 1000.0/(8.314472*300.0) # kJ/mol

# configuration

    # read initial configuration
    f = open(sys.argv[1], 'r')

    na = int(re.findall('\S+', f.readline())[0])
    xi = np.zeros((na), dtype=int)
    yi = np.zeros((na), dtype=int)
    zi = np.zeros((na), dtype=int)
    ti = np.zeros((na), dtype=int)

    ui = np.zeros((na), dtype=float)

    nx, ny, nz = map(int, re.findall('\S+', f.readline())[0:3])
    box = np.zeros((nx,ny,nz), dtype=int) - 1

    mv = []
    for i in range(na):
        ti[i], xi[i], yi[i], zi[i] = map(int, re.findall('\S+', f.readline())[0:4])
        box[xi[i], yi[i], zi[i]] = i
        if ti[i] == 1 or ti[i] == 2:
            mv.append(i)

    f.close()

    npm = len(mv)

    # initial energy
    for i in range(npm):
        j = mv[i] # select a movable particle (metal atoms)
        tti = ti[j]
        xt = xi[j]
        yt = yi[j]
        zt = zi[j]
        dj, du = upart_ini(j, box, ee, tti, xt, yt, zt,nx,ny,nz)
        ui[j] = sum(du)
    
    niter = 100000
    # nvt (uvt) simulation cycle 

    for it in range(niter):
        # move
        ia = 0
        for ir in range(npm):
            rm = rn.random()
            i = mv[rn.randint(0,npm-1)] # select a movable particle (metal atoms)
    
            upass = 0
            if rm < p_tran:
                mi = 0
                xt = (xi[i] + 2*(rn.randint(0,2) - 1)) % nx
                yt = (yi[i] + 2*(rn.randint(0,2) - 1)) % ny
                zt = zi[i]
                if box[xt,yt,zi[i]] != -1: # already occupied
                    upass = 1
                else:
                    tti = ti[i]
                    j = 0
            elif rm < p_mute:
                mi = 1
                tti = rn.randint(0,1) + 1 # Random type
                if tti == ti[i]:
                    upass = 1
                else:
                    xt = xi[i]
                    yt = yi[i]
                    zt = zi[i]
                    j = 0
            elif rm < p_swap:
                mi = 2
                j = mv[rn.randint(0,npm-1)] # select a movable particle
                if ti[i] == ti[j]:
                    upass = 1
                else:
                    tti = ti[j] # swap types
                    ttj = ti[i] # swap types
                    xt = xi[i]
                    yt = yi[i]
                    zt = zi[i]
                    xu = xi[j]
                    yu = yi[j]
                    zu = zi[j]
    
            if upass == 1:
                continue
    
            # energy difference calcualtion
            if mi == 0:
                dj, du = upart_t(i, box, ee, tti, xt, yt, zt,nx,ny,nz)
                ud = sum(du)

            if mi == 1 or mi == 2:
                dj, du = upart(i, box, ee, tti, xt, yt, zt,nx,ny,nz)
                ud = sum(du)
    
            if mi == 2:
                dj, du = upart(j, box, ee, ttj, xu, yu, zu,nx,ny,nz)
                ud = ud + sum(du)
    
            # condition
            if ud < 0.0: # accept
                accept(i, j, mi, dj, du, ui)
                ia = ia + 1
            elif np.exp(-beta*ud) > rn.random():
                accept(i, j, mi, dj, du, ui)
                ia = ia + 1

        if it % 10 == 0:
            cfgout(ti,xi,yi,zi,na,it)
            print it, np.sum(ui), float(ia)/float(npm), ia, npm

