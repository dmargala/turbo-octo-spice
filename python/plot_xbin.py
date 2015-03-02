#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.colors import LogNorm

import os
import sys

name = sys.argv[1]
prefix = '../output'

filename = os.path.join(prefix,'xibin_%s' % name)

if os.path.isfile(filename+'.npy'):
    print 'Reading %s.npy...' % filename
    xibin = np.load(filename+'.npy')
else:
    print 'Reading %s.txt...' % filename
    xibin = np.loadtxt(filename+'.txt')
    np.save(filename+'.npy', xibin)

Fi = xibin[:, 2]
Fj = xibin[:, 3]

print 'mean = ', np.mean(Fi), np.mean(Fj)
print 'Cov = ', np.cov(Fi, Fj)
print 'Corr = ', np.corrcoef(Fi, Fj)[0,1]

print 'transforming to Gaussian dist...'

def map2gauss(arr):
    random = np.random.normal(size=len(arr))
    arr[np.argsort(arr)] = random[np.argsort(random)]
map2gauss(Fi)
map2gauss(Fj)

print 'mean = ', np.mean(Fi), np.mean(Fj)
print 'Cov = ', np.cov(Fi, Fj)
print 'Corr = ', np.corrcoef(Fi, Fj)[0,1]

max_value = 5
num_bins = 50
bins = np.linspace(-max_value, max_value, num=num_bins+1)

fig = plt.figure(figsize=(8,6))
plt.hist(Fi, bins=bins, alpha=.4, color='red', label='$F_i$')
plt.hist(Fj, bins=bins, alpha=.4, color='blue', label='$F_j$')
plt.xlim(-max_value, max_value)
plt.ylabel('Counts')
plt.xlabel('$F_x$')
plt.grid()
plt.legend(loc=0)
plt.savefig(os.path.join(prefix,'bin_r%s_1d.png' % name))

fig = plt.figure(figsize=(8,6))
plt.hist2d(Fi, Fj, bins=bins, cmap='nipy_spectral_r', norm=LogNorm())
plt.title('$\Delta r \, \\approx \, %s $' % name)
plt.ylabel('$F_j$')
plt.xlabel('$F_i$')
plt.xlim(-max_value, max_value)
plt.ylim(-max_value, max_value)
plt.colorbar()
plt.savefig(os.path.join(prefix,'bin_r%s_2d.png' % name))
