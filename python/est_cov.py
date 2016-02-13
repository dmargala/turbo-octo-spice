#!/usr/bin/env python

import argparse
import glob
import numpy as np
import numpy.ma as ma

import os.path

from sklearn.covariance import empirical_covariance
from sklearn.covariance import MinCovDet
from sklearn.covariance import ledoit_wolf

def weighted_cov(X,W):
    """
    x,w shape is (num_variables, num_observations)
    """
    xx0 = X - ma.average(X, weights=W, axis=1)[:,np.newaxis]
    return ma.cov(xx0*np.sqrt(W))

def save_cov(fname, cov):
    num_entries = cov.shape[0]
    indices = np.arange(num_entries)
    xi, yi = np.meshgrid(indices,indices)
    tri_indices = np.tril_indices(num_entries)
    outdata = np.array([xi[tri_indices], yi[tri_indices], cov[tri_indices]])
    np.savetxt(fname, outdata.transpose(), fmt='%d %d %.18e')

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--name', type=str, default=None,
        help='filename glob pattern')
    parser.add_argument('--save', type=str, default='est.cov',
        help='output filename for cov estimate')
    args = parser.parse_args()

    filenames = glob.glob(args.name)
    if len(filenames) == 1:
        subsamples = np.load(filenames[0])
    else:
        first_xi = np.loadtxt(filenames[0])
        subsamples = np.empty((len(filenames), first_xi.shape[0], first_xi.shape[1]))
        for i,filename in enumerate(filenames):
            subsamples[i] = np.loadtxt(filename)
        np.save('cov_est_data', subsamples)

    num_samples, num_bins, num_dim = subsamples.shape

    xis = subsamples[:,:,1]
    wgts = subsamples[:,:,2]

    print xis.shape

    masked_wgts = ma.MaskedArray(wgts, mask=(wgts == 0))
    masked_xis = ma.MaskedArray(xis, mask=(wgts == 0))
    cov = weighted_cov(masked_xis.T, masked_wgts.T) / 6.4798e8

    print cov.shape

    # check symmetric
    if not np.allclose(cov.transpose(), cov):
        raise RuntimeError('Covariance estimate is not symmetric')

    # check if pos def
    try:
        np.linalg.cholesky(cov)
    except np.linalg.LinAlgError:
        print 'Covariance is not positive definite, attempting to fix...'
        w,v = np.linalg.eig(cov.data)
        w[w < 0] = 1e-14
        cov = v.dot(np.diagflat(w).dot(np.linalg.inv(v)))
        if not np.isreal(cov).all():
            raise ValueError('Covariance has complex elements!')

    print 'Saving covariance matrix...'

    print '(sign, logdet) : ', np.linalg.slogdet(cov)

    save_cov(args.save, cov)


if __name__ == '__main__':
    main()
