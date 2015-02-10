#!/usr/bin/env python
"""
"""

import argparse

import h5py
import healpy as hp
import numpy as np

z0 = 2.5 
scale = 3802.63
qso_z = 2.7

redshifts = [2.25569, 2.2614, 2.26712, 2.27285, 2.2786, 2.28436, 2.29014, 2.29593, \
2.30173, 2.30755, 2.31339, 2.31924, 2.3251, 2.33098, 2.33687, \
2.34278, 2.34871, 2.35464, 2.3606, 2.36656, 2.37255, 2.37855, \
2.38456, 2.39059, 2.39663, 2.40269, 2.40877, 2.41486, 2.42096, \
2.42708, 2.43322, 2.43937, 2.44554, 2.45173, 2.45793, 2.46414, \
2.47037, 2.47662, 2.48289, 2.48917, 2.49546, 2.50177, 2.5081, \
2.51445, 2.52081, 2.52719, 2.53358, 2.54, 2.54642, 2.55287, 2.55933, \
2.56581, 2.5723, 2.57882, 2.58534, 2.59189, 2.59845, 2.60503, \
2.61163, 2.61824, 2.62488, 2.63153, 2.63819, 2.64488]

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
        help='more verbose output')
    parser.add_argument('--input', type=str, default=None,
        help='input file name')
    parser.add_argument('--nx', type=int, default=None,
        help='pts per axis')
    parser.add_argument('--spacing', type=float, default=5,
        help='grid point spacing in Mpc/h')
    parser.add_argument('--output', type=str, default=None,
        help='output file name')
    args = parser.parse_args()

    dtheta = args.spacing/scale

    data = np.loadtxt(args.input)
    reshaped = data.reshape(args.nx, args.nx, args.nx, 5)

    outfile = h5py.File(args.output,'w')
    delta_field = outfile.create_group('delta_field')
    for ix in range(args.nx):
        for iy in range(args.nx):
            los_grp = delta_field.create_group('%d-%d' % (ix, iy))
            los_grp.attrs['ra'] = ix*dtheta*180.0/np.pi
            los_grp.attrs['dec'] = iy*dtheta*180.0/np.pi
            los_grp.attrs['z'] = qso_z
            los_grp.create_dataset('absorber_z', data=redshifts, dtype='f4')
            los_grp.create_dataset('absorber_delta', data=reshaped[ix,iy,:,3], dtype='f4')
            los_grp.create_dataset('absorber_ivar', data=reshaped[ix,iy,:,4], dtype='f4')






if __name__ == '__main__':
    main()