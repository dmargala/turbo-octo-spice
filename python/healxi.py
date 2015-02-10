#!/usr/bin/env python
"""
"""

import argparse

import h5py
import healpy as hp
import numpy as np

import matplotlib.pyplot as plt


def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--verbose', action='store_true',
        help='more verbose output')
    parser.add_argument('--order', type=int, default=5,
        help='order of healpix hierachical tesselation')
    parser.add_argument('--input', type=str, default=None,
        help='input file name')
    args = parser.parse_args()

    nside = 2**args.order
    npix = hp.nside2npix(nside)
    m = np.arange(npix)

    print 'HEALPix map info:'
    print '\tn pixels:', npix
    print '\tresolution:', hp.nside2resol(nside)
    print '\tpixel area:', hp.nside2pixarea(nside)
    print '\tmax distance from center:', hp.max_pixrad(nside)

    print 'xi binning info:'
    transerve_comoving_scale = 3820.43 # Mpc, at z = 2.1
    maxang = 200/transerve_comoving_scale
    print '\tmax ang scale:', maxang

    print 'test stuff:'

    test_theta = 5./6.*np.pi/2
    test_phi = 7./6.*np.pi
    testvec = hp.ang2vec(test_theta, test_phi)

    print '\ttheta, phi:', test_theta, test_phi
    print '\t3d vec:', testvec
    print '\tpixel index:', hp.ang2pix(nside, test_theta, test_phi)
    print '\tclosest neighbors:', hp.pixelfunc.get_all_neighbours(nside, test_theta, test_phi)
    print '\tn neighbors within max ang scale:', len(hp.query_disc(nside, testvec, maxang, inclusive=True))

    include_pix = hp.query_disc(nside, hp.ang2vec(5./6.*np.pi/2, 7./6.*np.pi), 2*maxang, inclusive=True)

    print '\tn pixels for test:', len(include_pix)

    qso_pos = np.loadtxt(args.input)
    print 'Shape of input data: ', qso_pos.shape

    m = {}

    qso_pos[:,0] = np.deg2rad(qso_pos[:,0])
    qso_pos[:,1] = np.deg2rad(90-qso_pos[:,1])

    class Quasar():
        def __init__(self, theta, phi, z):
            self.theta = theta
            self.phi = phi
            self.z = z
            # fake this for now
            self.wave = np.log10((1+z)*(1100+np.arange(10)))

    quasars = []

    for i, qso in enumerate(qso_pos):
        ra, dec, z = qso
        phi = ra
        theta = dec
        q = Quasar(theta, phi, z)
        map_index = hp.ang2pix(nside, theta, phi)

        if map_index not in include_pix:
            continue
        #q.map_index = map_index
        if map_index in m:
            m[map_index].append(i)
        else:
            m[map_index] = [i]
        q.i = i
        quasars.append(q)

    print 'Number of pixels with at least one quasar: %d' % (len(m.keys()))
    print 'Average number of quasars per pixel: %f' % (float(len(quasars))/len(m.keys()))

    nearby_pixels_sum = 0
    nearby_quasars_sum = 0
    for q_index, q in enumerate(quasars):
        nearby_pixels = hp.query_disc(nside, hp.ang2vec(q.theta, q.phi), maxang, inclusive=True)
        nearby_pixels_sum += len(nearby_pixels)
        for nearby_pixel in nearby_pixels:
            if nearby_pixel in m:
                nearby_quasars = m[nearby_pixel]
                nearby_quasars_sum += len(nearby_quasars)

                p_phi = qso_pos[nearby_quasars,0]
                p_theta = qso_pos[nearby_quasars,1]

                #p1.sth*p2.sth*(p1.cph*p2.cph + p1.sph*p2.sph) + p1.cth*p2.cth;
                angdist = np.sin(q.theta)*np.sin(p_theta)*(np.cos(q.phi)*np.cos(p_phi) + np.sin(q.phi)*np.sin(p_phi)) + np.cos(q.theta)*np.cos(p_theta)

                for p_i, p_index in enumerate(nearby_quasars):
                    theta, phi, z = qso_pos[p_index]
                    p = Quasar(theta, phi, z)
                    loglamdist = abs(np.subtract.outer(q.wave, p.wave))
                    p_angdist = angdist[p_i]
                #losdist = abs(np.log(q_wave/p.wave))
    print 'Average number of pixels within 200 Mpc of a quasar: %f' % (float(nearby_pixels_sum)/len(quasars))
    print 'Average number of quasars within nearby pixels of a quasar: %f' % (float(nearby_quasars_sum)/len(quasars))

    mask = np.ones(npix).astype(np.bool)
    mask[m.keys()] = False

    mm = hp.ma(np.zeros(npix))
    mm.mask = mask
    mm[~mm.mask] = [len(m[i]) for i in m.keys()]

    hp.mollview(mm.filled(), rot=(270,0,0), flip='geo')

    #plt.show()



if __name__ == '__main__':
    main()