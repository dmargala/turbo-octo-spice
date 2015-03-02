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
    filename = os.path.join(prefix,'healxi-%s' % name)
    index, r, mu, z, didj, di, dj, wgt = np.loadtxt(filename+'.txt', unpack=True)
    err = 1.0/np.sqrt(wgt)
    return r, didj, err

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

r, xi1024, err1024 = loadXi('1024-5')
r, xi512, err512 = loadXi('512-10')
r, xi512b, err512b = loadXi('512-10.2')

fig = plt.figure(figsize=(8,6))

plt.plot(r_fid, r_fid**2*(xi_fid/4.0), 'r', label='PlanckWPBestFitLCDM')

plt.plot(r, r**2*(xi512), 'b', label='cosmo grf (512/10)')
plt.errorbar(r, r**2*(xi512), yerr=r**2*err512, color='b')

plt.plot(r, r**2*(2*xi512b), 'b--', label='cosmo grf (512/10)')

plt.plot(r, r**2*(xi1024), 'g', label='cosmo grf (1024/5)')
plt.errorbar(r, r**2*(xi1024), yerr=r**2*err1024, color='g')

plt.legend()

plt.xlim(0, 200)
plt.ylabel('$r^2 \\xi(|r|)$')
plt.xlabel('Comoving Separation $|r|$ (Mpc/h)')
plt.grid(ls='-', alpha=.3)
plt.savefig(outfilename+'.png')

