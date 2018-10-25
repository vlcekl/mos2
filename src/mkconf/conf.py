#!/usr/bin/python
#
# File name:   plot.py
# Date:        2007/01/15 15:05
# Author:      Lukas Vlcek
#
# Description: 
#

import matplotlib
#matplotlib.use('PS')
from pylab import *
from mpl_toolkits.mplot3d import axes3d, Axes3D
from matplotlib.font_manager import *
from matplotlib.mlab import load
import scipy.io as sio
from scipy.fftpack import ifftn
import re

def findmax(img, ai, aj):
    ima = img.shape[0]
    d = img[ai,aj]
    mxi = ai
    mxj = aj
    for k in range(-1,2):
        ki = ai + k
        if ki < 0 or ki > ima-1:
            continue
        for l in range(-1,2):
            lj = aj + l
            if lj < 0 or lj > ima-1:
                continue
            if img[ki,lj] >= d:
                d = img[ki,lj]
                mxi = ki
                mxj = lj
    return (mxi, mxj)


mm = sio.loadmat('good_image.mat')

im = np.array(mm['tophat_image'])
at = np.array(mm['A'])

ima = plt.imshow(im[:, :], origin='lower', cmap=cm.binary)

draw()
savefig('im_g')

# create grid
# analyze densities
# assign pixels to atoms/micro-regions
nm = im.shape[0]
rg = np.zeros((nm, nm), dtype=float) # region assignement
rm = []

kr = 0
for i in range(nm):
    for j in range(nm):
        ai = i
        aj = j

        if rg[ai,aj] != 0:
            continue
        pth = []

        pth.append([ai, aj])
        (mi, mj) = findmax(im, ai, aj)

        while mi != ai or mj != aj:
            if rg[mi,mj] != 0:  # if the max neighbor belongs to a region, add the new path to the region
                for pt in pth:
                    rg[pt[0], pt[1]] = rg[mi,mj]
                pth = []
                break
            else:               # else add the max neighbor to the current path and find a new max
                pth.append([mi, mj])
                ai = mi
                aj = mj
                (mi, mj) = findmax(im, ai, aj)


        if len(pth) > 0:    
            kr = kr + 1
            rm.append([mi, mj])
            for pt in pth:
                rg[pt[0], pt[1]] = kr

mkr = kr

rb = np.zeros((nm, nm), dtype=float) # region assignement
for i in range(nm):
    for j in range(nm):
        for k in range(-1,2):
            ki = i + k
            if ki < 0 or ki > nm-1:
                continue
            for l in range(-1,2):
                lj = j + l
                if lj < 0 or lj > nm-1:
                    continue
                if rg[ki,lj] != rg[i,j]:
                    rb[i,j] = 1
                    break

imb = plt.imshow(rb[:, :], origin='lower', cmap=cm.binary)

draw()
savefig('im_b')

# COMs and intensities of regions
krs = np.zeros((mkr), dtype=float)
rrx = np.zeros((mkr), dtype=float)
rry = np.zeros((mkr), dtype=float)

print krs.shape, im.shape, rrx.shape, rb.shape

for i in range(nm):
    for j in range(nm):
        k = rg[i,j]-1

        krs[k] = krs[k] + im[i,j]
        rrx[k] = rrx[k] + im[i,j]*(float(j)+0.5)
        rry[k] = rry[k] + im[i,j]*(float(i)+0.5)

rrx[:] = rrx[:]/krs[:]
rry[:] = rry[:]/krs[:]

print krs
print rrx
print rry

# region intensity spatial assignement
ri = np.zeros((nm, nm), dtype=float)
for i in range(nm):
    for j in range(nm):
        ri[i,j] = krs[rg[i,j]-1]

imi = plt.imshow(ri[:, :], origin='lower', cmap=cm.binary)

draw()
savefig('im_i')

# intensity filter to show only significant COMs
flsx = []
flsy = []
mass = []

for k in range(mkr):
    #if krs[k] > 1.0:
    if krs[k] > 1.0:
        if rrx[k] > 5 and rrx[k] < 508 and rry[k] > 5 and rry[k] < 508:
            flsx.append(rrx[k])
            flsy.append(rry[k])
            mass.append(krs[k])

scatter(flsx[:], flsy[:], c='r',s=5)

draw()
savefig('im_s')

lmax = im.shape[0]

figure(1, figsize =(8.3, 6.0), dpi=600)

plt.imshow(im[:, :], origin='lower', cmap=cm.binary)

mx = np.zeros((lmax), dtype=float)
my = np.zeros((lmax), dtype=float)
ny = np.zeros((lmax), dtype=float)
for i in range(lmax):
    my[i] = float(i)
    #my[i] = float(i)
    mx[i] = -float(i)*np.sin(2.7*np.pi/180.0)
    ny[i] = -float(i)*np.sin(0.63*np.pi/180.0)
#    ny[i] = 0.0

for i in range(54):
    #fi = float(i - 27)*9.99
    #fj = float(i - 27)*9.85
    #plot(my[:], mx[:]+268.6+fi+4.5,'b',linewidth=0.1)
    #plot(ny[:]+259.5+fj-1.0, my[:],'r',linewidth=0.1)
    fi = float(i - 27)*10.012
    fj = float(i - 27)*9.862
    plot(my[:], mx[:]+268.6+fi+4.5+0.05,'b',linewidth=0.3)
    plot(ny[:]+259.5+fj-1.0-1.2, my[:],'r',linewidth=0.3)

scatter(flsx[:], flsy[:], c='r', marker='s', s=30.0)

axis([0,512,0,512])

draw()
savefig('all')

# analyze particle positions within the lattice

# determine lattice unit vectors
ax = 1.0
ay = (-np.sin(2.7*np.pi/180.0))

bx = (-np.sin(0.63*np.pi/180.0))
by = 1.0

aa = (ax*ax + ay*ay)**0.5
a0x = ax/aa
a0y = ay/aa

bb = (bx*bx + by*by)**0.5
b0x = bx/bb
b0y = by/bb

#da = 9.85*a0x
da = 9.862*a0x
#db = 9.99*b0y
db = 10.012*b0y

# angle between unit vectors

cosa = a0x*b0x + a0y*b0y
sina = (1.0 - cosa**2)**0.5
print 'kozina', cosa, sina

clsx = []
clsy = []
mlsx = []
mlxy = []

dac = 268.35 - int(268.35/9.85)*9.85 + 0.05
dbc = 283.09 - int(283.09/9.99)*9.99 - 1.2

# determine COM position in lattice coordinates

latx = np.zeros((100, 100), dtype=float) # lattice sites
laty = np.zeros((100, 100), dtype=float) # lattice sites
lata = np.zeros((100, 100), dtype=float) # lattice sites
latb = np.zeros((100, 100), dtype=float) # lattice sites
latm = np.zeros((100, 100), dtype=float) # lattice sites
#lati = np.zeros((100, 100), dtype=float) # lattice sites

iamax = 0
ibmax = 0
iamin = 0
ibmin = 0
for i in range(len(flsx)):
    xx = flsx[i] - dac
    yy = flsy[i] - dbc
    #r = (flsx[i]**2 + flsy[i]**2)**0.5
    #d = flsx[i]*a0x + flsy[i]*a0y
    r = (xx**2 + yy**2)**0.5
    d = xx*a0x + yy*a0y
    y = (r*r - d*d)**0.5

    br = y/sina
    if cosa < 0:
        ar = d + (br*br - y*y)**0.5
    else:
        ar = d - (br*br - y*y)**0.5

    iar = int(ar/da)
    ibr = int(br/db)
    dar = ar - float(iar)*da
    dbr = br - float(ibr)*db

    if iar > iamax:
        iamax = iar
    if ibr > ibmax:
        ibmax = ibr
    if iar < iamin:
        iamin = iar
    if ibr < ibmin:
        ibmin = ibr

    print i, flsx[i], flsy[i], ar, br, iar, ibr, iar*da, ibr*db, dar, dbr

    lata[iar, ibr] = lata[iar, ibr] + dar*mass[i]
    latb[iar, ibr] = latb[iar, ibr] + dbr*mass[i]
    #latx[iar, ibr] = latx[iar, ibr] + flsx[i]*mass[i]
    #laty[iar, ibr] = laty[iar, ibr] + flsy[i]*mass[i]
    #latx[iar, ibr] = (float(iar)+0.5)*a0x*9.862 + (float(ibr)+0.5)*bx*10.012 + dac
    #laty[iar, ibr] = (float(iar)+0.5)*a0y*9.862 + (float(ibr)+0.5)*by*10.012 + dbc
    #latx[iar, ibr] = (float(iar)+0.5)*a0x*9.862 + (float(ibr)+0.5)*b0x*9.862 + dac
    #laty[iar, ibr] = (float(iar)+0.5)*a0y*10.012 + (float(ibr)+0.5)*b0y*10.012 + dbc
    latx[iar, ibr] = (float(iar)+0.5)*a0x*da*1.0028 + (float(ibr)+0.5)*b0x*db + dac - 0.5
    laty[iar, ibr] = (float(iar)+0.5)*a0y*da + (float(ibr)+0.5)*b0y*db + dbc
    latm[iar, ibr] = latm[iar, ibr] + mass[i]

#    xx = ar*a0x + br*b0x
#    yy = ar*a0y + br*b0y
#
#    print i, flsx[i], flsy[i], ar, br, r, (ar*ar + br*br)**0.5, (xx*xx+yy*yy)**0.5

print da, db

lati = []
for iar in range(iamin, iamax+1):
    for ibr in range(ibmin, ibmax+1):
        latx[iar, ibr] = (float(iar)+0.5)*a0x*da*1.0028 + (float(ibr)+0.5)*b0x*db + dac - 0.5 - 0.0
        laty[iar, ibr] = (float(iar)+0.5)*a0y*da + (float(ibr)+0.5)*b0y*db + dbc + 1.0
        #latx[iar, ibr] = (float(iar)+0.5)*a0x*da + (float(ibr)+0.5)*b0x*db + dac
        #laty[iar, ibr] = (float(iar)+0.5)*a0y*da + (float(ibr)+0.5)*b0y*db + dbc
        if latm[iar, ibr] > 0.0:
            lata[iar, ibr] = lata[iar, ibr]/latm[iar, ibr]
            latb[iar, ibr] = latb[iar, ibr]/latm[iar, ibr]
#            latx[iar, ibr] = latx[iar, ibr]/latm[iar, ibr]
#            laty[iar, ibr] = laty[iar, ibr]/latm[iar, ibr]
            lati.append([float(iar), float(ibr), latm[iar,ibr], latx[iar,ibr], laty[iar,ibr]])

lati = np.array(lati)
figure(2, figsize =(6.0, 6.0), dpi=600)

#for i in range(-4,5):
#    for j in range(-4,5):
#for i in range(3,4):
#    for j in range(-1,0):
#        ia = int(latx[20,20])
#        ib = int(laty[20,20])
#        im[ib+j, ia+i] = 1.2
#        ia = int(latx[30,30])
#        ib = int(laty[30,30])
#        im[ib+j, ia+i] = 1.2
#        ia = int(latx[40,40])
#        ib = int(laty[40,40])
#        im[ib+j, ia+i] = 1.2
#        ia = int(latx[20,40])
#        ib = int(laty[20,40])
#        im[ib+j, ia+i] = 1.2
#
#for iar in range(iamin, iamax+1):
#    for ibr in range(ibmin, ibmax+1):
#        ii = int(latx[iar,ibr])
#        jj = int(laty[iar,ibr])
#        if ii >= 0 and ii < 512 and jj >= 0 and jj < 512:
#            im[jj,ii] = 1.2

plt.imshow(im[:, :], interpolation='none', origin='lower', cmap=cm.binary)
#scatter(lati[:,0], lati[:,1], c='r', marker=(4,0), s=35.0)
#scatter(lati[:,3], lati[:,4], c='r', marker='s', s=35.0)
scatter(lati[:,3], lati[:,4], c='r', marker='s', s=0.1)
#axis([iamin,iamax,ibmin,ibmax])
for i in range(54):
    fi = float(i - 27)*10.012
    fj = float(i - 27)*9.862
    plot(my[:], mx[:]+268.6+fi+4.5+0.05,'b',linewidth=0.3)
    plot(ny[:]+259.5+fj-1.0-1.2, my[:],'r',linewidth=0.3)
axis([0,512,0,512])

draw()
savefig('latt')

figure(4, figsize =(6.0, 6.0), dpi=600)
scatter(lati[:,0], lati[:,1], c='r', marker='s', s=35)
axis([iamin,iamax,ibmin,ibmax])
draw()
savefig('chess')

r1 = (da + db)/4.0
r2 = r1*2.0**0.5
ns = int(r1)+1
#r1 = 0.15*r1
ns = 4
r1 = 4
r2 = 5
ns = 4

sbx = np.zeros((100, 100), dtype=float) # lattice sites

s_min = np.zeros((2*ns+1, 2*ns+1), dtype=float) # lattice sites
s_sd1 = np.zeros((2*ns+1, 2*ns+1), dtype=float) # lattice sites
s_sd2 = np.zeros((2*ns+1, 2*ns+1), dtype=float) # lattice sites
s_sd3 = np.zeros((2*ns+1, 2*ns+1), dtype=float) # lattice sites
s_sd4 = np.zeros((2*ns+1, 2*ns+1), dtype=float) # lattice sites

# masks
m_min = np.zeros((2*ns+1, 2*ns+1), dtype=float)
m_sd1 = np.zeros((2*ns+1, 2*ns+1), dtype=float)
m_sd2 = np.zeros((2*ns+1, 2*ns+1), dtype=float)
m_sd3 = np.zeros((2*ns+1, 2*ns+1), dtype=float)
m_sd4 = np.zeros((2*ns+1, 2*ns+1), dtype=float)

m_min[2:7,2:7] = 1.0
m_sd1[0,1:8] = 1.0
m_sd2[1:8,0] = 1.0
m_sd3[8,1:8] = 1.0
m_sd4[1:8,8] = 1.0

for iar in range(iamin, iamax+1):
    for ibr in range(ibmin, ibmax+1):
        ia = int(latx[iar, ibr])
        ib = int(laty[iar, ibr])
        s_min[:,:] = 0.0
        s_sd1[:,:] = 0.0
        s_sd2[:,:] = 0.0
        s_sd3[:,:] = 0.0
        s_sd4[:,:] = 0.0
        for i in range(-ns,ns+1):
            for j in range(-ns,ns+1):
        #for i in range(-ns,ns):
        #    for j in range(-ns,ns):
                ii = ib+j
                jj = ia+i
                if ii >= 0 and ii < 512 and jj >= 0 and jj < 512:
                    sbx[ibr, iar] = sbx[ibr, iar] + im[ii, jj]
                    s_min[i+ns,j+ns] = s_min[i+ns,j+ns] + im[ii, jj]*m_min[i+ns,j+ns]
                    s_sd1[i+ns,j+ns] = s_sd1[i+ns,j+ns] + im[ii, jj]*m_sd1[i+ns,j+ns]
                    s_sd2[i+ns,j+ns] = s_sd2[i+ns,j+ns] + im[ii, jj]*m_sd2[i+ns,j+ns]
                    s_sd3[i+ns,j+ns] = s_sd3[i+ns,j+ns] + im[ii, jj]*m_sd3[i+ns,j+ns]
                    s_sd4[i+ns,j+ns] = s_sd4[i+ns,j+ns] + im[ii, jj]*m_sd4[i+ns,j+ns]

#        if np.sum(s_min) < np.sum(s_sd1):#*25./18.:
#            sbx[ibr, iar] = 0.0
#        if np.sum(s_min) < np.sum(s_sd2):#*25./18.:
#            sbx[ibr, iar] = 0.0
#        if np.sum(s_min) < np.sum(s_sd3):#*25./18.:
#            sbx[ibr, iar] = 0.0
#        if np.sum(s_min) < np.sum(s_sd4):#*25./18.:
#            sbx[ibr, iar] = 0.0
        if np.sum(s_min) < np.sum(s_sd4)+np.sum(s_sd3)+np.sum(s_sd2)+np.sum(s_sd1):#*25./18.:
            sbx[ibr, iar] = 0.0

        if (ibr+iar)%2 == 0:
            sbx[ibr, iar] = -1.0*sbx[ibr, iar]


#                rr = float(ii*ii + jj*jj)**0.5
#                if irx+ii<0 or irx+ii>511 or iry+jj<0 or iry+jj>511:
#                    continue
#                print rx, ry,irx, iry, irx+ii, iry+jj
#                if rr < r2:
#                    if ri < r1 and rj < r1:
#                        if rr < r1: # small circle
#                            sr1[iar, ibr] = sr1[iar, ibr] + im[irx+ii, iry+jj]
#                        else: #box
#                            sbx[iar, ibr] = sbx[iar, ibr] + im[irx+ii, iry+jj]
#                    else: # large circle
#                        sr2[iar, ibr] = sr2[iar, ibr] + im[irx+ii, iry+jj]
#for i in range(-4,5):
#    for j in range(-4,5):
#for i in range(3,4):
#    for j in range(-1,0):
#        ia = int(latx[20,20])
#        ib = int(laty[20,20])
#        im[ib+j, ia+i] = 1.2

figure(3, figsize =(6.0, 6.0), dpi=600)

smin = np.amin(sbx[iamin:iamax+1, ibmin:ibmax+1])
smax = np.amax(sbx[iamin:iamax+1, ibmin:ibmax+1])
if smax > -smin:
    sbx[iamax,ibmax] = -smax
else:
    sbx[iamax,ibmax] = -smin

smin = np.amin(sbx[iamin:iamax+1, ibmin:ibmax+1])
smax = np.amax(sbx[iamin:iamax+1, ibmin:ibmax+1])

print 'minimax', smin, smax

plt.imshow(sbx[iamin:iamax+1, ibmin:ibmax+1], interpolation='none', origin='lower', cmap=cm.bwr)
axis([iamin,iamax,ibmin,ibmax])

draw()
savefig('lattint')

figure(5, figsize =(6.0, 6.0), dpi=600)

sb = np.zeros((100, 100), dtype=float) # lattice sites
sn = np.zeros((100, 100), dtype=float) # lattice sites

for iar in range(iamin, iamax+1):
    for ibr in range(ibmin, ibmax+1):
        cnt = 0
        if sbx[iar, ibr] < 0.0:
            if sbx[iar+1,ibr+1]<0. and sbx[iar,ibr+1]==0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar+1,ibr-1]<0. and sbx[iar,ibr-1]==0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-1,ibr-1]<0. and sbx[iar,ibr-1]==0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-1,ibr+1]<0. and sbx[iar,ibr+1]==0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar+2,ibr]<0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-2,ibr]<0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr+2]<0. and sbx[iar,ibr+1]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr-2]<0. and sbx[iar,ibr-1]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr-1] > 0. and sbx[iar,ibr+1] > 0.:
                cnt = 0
            if sbx[iar-1,ibr] > 0. and sbx[iar+1,ibr] > 0.:
                cnt = 0
        elif sbx[iar, ibr] > 0.0:
            if sbx[iar+1,ibr+1]>0. and sbx[iar,ibr+1]==0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar+1,ibr-1]>0. and sbx[iar,ibr-1]==0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-1,ibr-1]>0. and sbx[iar,ibr-1]==0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-1,ibr+1]>0. and sbx[iar,ibr+1]==0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar+2,ibr]>0. and sbx[iar+1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar-2,ibr]>0. and sbx[iar-1,ibr]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr+2]>0. and sbx[iar,ibr+1]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr-2]>0. and sbx[iar,ibr-1]==0.:
                cnt = cnt + 1
            if sbx[iar,ibr-1] < 0. and sbx[iar,ibr+1] < 0.:
                cnt = 0
            if sbx[iar-1,ibr] < 0. and sbx[iar+1,ibr] < 0.:
                cnt = 0
        else:
            cnt = -1

        if cnt > 0:
            sb[iar,ibr] = sbx[iar,ibr]
        elif cnt == 0:
            sn[iar,ibr] = sbx[iar,ibr]

smin = np.amin(sb[iamin:iamax+1, ibmin:ibmax+1])
smax = np.amax(sb[iamin:iamax+1, ibmin:ibmax+1])
if smax > -smin:
    sb[iamax,ibmax] = -smax
else:
    sb[iamax,ibmax] = -smin

smin = np.amin(sb[iamin:iamax+1, ibmin:ibmax+1])
smax = np.amax(sb[iamin:iamax+1, ibmin:ibmax+1])

print 'minimax', smin, smax

sb = np.zeros((100, 100), dtype=float) # lattice sites

dat = []
f = open('iii', 'r')
for line in f:
    sarr = re.findall('\S+', line)
    ix = int(sarr[0])
    iy = int(sarr[1])
    iii = float(sarr[2])
    if ix < 100 and iy < 100:
        sb[iy,ix] = -iii
f.close()

sb[iamax,ibmax] = -smax
sb[iamin,ibmin] = smax

plt.imshow(sb[iamin:iamax+1, ibmin:ibmax+1], interpolation='none', origin='lower', cmap=cm.bwr)
axis([iamin,iamax,ibmin,ibmax])

for i in range(iamin, iamax+1):
    for j in range(ibmin, ibmax+1):
        print sb[i,j]
#        if np.abs(sb[i,j]) > 1.e-5:
#            print 'O', j*2, i*2, 1

draw()
savefig('latgood_sim')

figure(6, figsize =(6.0, 6.0), dpi=600)

#smin = np.amin(sn[iamin:iamax+1, ibmin:ibmax+1])
#smax = np.amax(sn[iamin:iamax+1, ibmin:ibmax+1])
#if smax > -smin:
#    sn[iamax,ibmax] = -smax
#else:
#    sn[iamax,ibmax] = -smin
#
#smin = np.amin(sn[iamin:iamax+1, ibmin:ibmax+1])
#smax = np.amax(sn[iamin:iamax+1, ibmin:ibmax+1])

sn = np.abs(sn)

plt.imshow(sn[iamin:iamax+1, ibmin:ibmax+1], interpolation='none', origin='lower', cmap=cm.binary)
axis([iamin,iamax,ibmin,ibmax])

draw()
savefig('latbad')

figure(7, figsize =(6.0, 6.0), dpi=600)
sint = np.zeros((9), dtype=float)
mint = np.zeros((9), dtype=float)

for iar in range(iamin, iamax+1):
    for ibr in range(ibmin, ibmax+1):
        if sb[iar,ibr] == 0.0:
            continue
        cnt = 0
#        if sbx[iar,ibr+1] != 0.0:
#            cnt = cnt + 1
#        if sbx[iar+1,ibr] != 0.0:
#            cnt = cnt + 1
#        if sbx[iar-1,ibr] != 0.0:
#            cnt = cnt + 1
#        if sbx[iar,ibr-1] != 0.0:
#            cnt = cnt + 1
        if sb[iar+1,ibr+1] != 0.0:
            cnt = cnt + 1
        if sb[iar+1,ibr-1] != 0.0:
            cnt = cnt + 1
        if sb[iar-1,ibr-1] != 0.0:
            cnt = cnt + 1
        if sb[iar-1,ibr+1] != 0.0:
            cnt = cnt + 1

        mint[cnt] = mint[cnt] + 1
        sint[cnt] = sint[cnt] + np.abs(sb[iar,ibr])

for i in range(9):
    if mint[i] > 0:
        sint[i] = sint[i]/mint[i]

#plot(range(5), sint[:5], 'r',linewidth=0.3)
#plot(range(5), mint[:5]/np.sum(mint[:5]), 'b',linewidth=0.3)

#bar(range(5), sint[:5], color='r',linewidth=0.3)
bar(np.arange(5)-0.4, mint[:5]/np.sum(mint[:5]), color='b')
axis([-0.5, 4.5, 0.0, 0.4])

draw()
savefig('nbr')
figure(8, figsize =(6.0, 6.0), dpi=600)
bar(np.arange(5)-0.4, sint[:5], color='r')
axis([-0.5, 4.5, 0.0, 2.4])
draw()
savefig('nint')
