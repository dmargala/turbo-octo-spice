#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

import os
import sys

prefix = 'output'


if len(sys.argv) > 1:
    outfilename = os.path.join(prefix, sys.argv[1])
else:
    outfilename = os.path.join(prefix, 'corr1d')

def loadXi(name):
    filename = os.path.join(prefix, name)
    index, r, mu, z, didj, di, dj, wgt = np.loadtxt(filename+'.txt', unpack=True)
    err = 1.0/np.sqrt(wgt)
    return r, didj, err

def loadXiSimple(name):
    filename = os.path.join(prefix, name)
    index, didj = np.loadtxt(filename+'.txt', unpack=True)
    return index, didj


def C_l(ell, beta):
    if ell == 0:
        return 1 + 2.0/3.0*beta + 1.0/5.0*beta**2
    elif ell == 2:
        return 4.0/3.0*beta + 4.0/7.0*beta**2
    elif ell == 4:
        return 8.0/35.0*beta**3
    else:
        return 0

def loadMultipole(ell):
    model_filename = '../baofit/models/PlanckWPBestFitLCDM.%d.dat' % ell
    return np.loadtxt(model_filename, unpack=True)

r_fid, xi0_fid = loadMultipole(0)
r_fid, xi2_fid = loadMultipole(2)
r_fid, xi4_fid = loadMultipole(4)
beta = 0
xi_fid = C_l(0, beta)*xi0_fid + C_l(2, beta)*xi2_fid + C_l(4, beta)*xi4_fid

r, xi512, err512 = loadXi('healxi-512-10')
rtest, xitest, errtest = loadXi('healxi-256-2')
gpuindex, gpuxi = loadXiSimple('gpulos-xi-256-2')
rgpu = 4*(gpuindex+.5)

fig = plt.figure(figsize=(8, 6))

plt.plot(r_fid, r_fid**2*(xi_fid), 'r', label='PlanckWPBestFitLCDM')

plt.plot(rtest, rtest**2*(xitest), 'b', label='heal xi (256-2)')
# plt.errorbar(r, r**2*(xi512test), yerr=r**2*err512test, color='b')

plt.plot(rgpu, rgpu**2*(gpuxi), 'b', label='gpu xi los (256-2)')
# plt.errorbar(r, r**2*(xi512), yerr=r**2*err512, color='b')

plt.legend()

plt.xlim(0, 200)
plt.ylim(-5,10)
plt.ylabel('$r^2 \\xi(|r|)$')
plt.xlabel('Comoving Separation $|r|$ (Mpc/h)')
plt.grid(ls='-', alpha=.3)
plt.savefig(outfilename+'.png')

