#!/usr/bin/env python

import argparse
import glob
import numpy as np

import os.path

from sklearn.covariance import empirical_covariance
from sklearn.covariance import MinCovDet
from sklearn.covariance import ledoit_wolf

def main():
    # parse command-line arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--name', type=str, default=None,
        help='filename glob pattern')
    parser.add_argument('--cache', action='store_true',
        help='read input data from cache')
    parser.add_argument('--save', type=str, default='est.cov',
        help='output filename for cov estimate')
    args = parser.parse_args()

    if args.cache and os.path.isfile(args.name):
        all_xis = np.load(args.name)
    else:
        filenames = glob.glob(args.name)

        num_xis = len(filenames)

        first_xi = np.loadtxt(filenames[0])
        print first_xi.shape

        all_xis = np.empty((num_xis, first_xi.shape[0], first_xi.shape[1]))

        for i,filename in enumerate(filenames):
            all_xis[i] = np.loadtxt(filename)

        np.save('cov_est_data', all_xis)
    num_samples, num_bins, num_dim = all_xis.shape

    mean_xi = np.average(all_xis[:,:,1], weights=all_xis[:,:,2], axis=0)
    sum_weights = np.sum(all_xis[:,:,2], axis=0)

    cov = 1.0/sum_weights

    with open(args.save, 'w') as outfile:
        for i, value in enumerate(cov):
            outfile.write('{} {} {}\n'.format(i,i,value))

#    my_cov_est = np.average((all_xis[:,:,1] - mean_xi)**2, weights=all_xis[:,:,2], axis=0)

    #X = all_xis[:,:,1]
    #cov = ledoit_wolf(X, assume_centered=True)[0]

    #mcd = MinCovDet(assume_centered=True).fit(X)
    #print mcd

    print mean_xi
    #print cov


    #print np.linalg.det(cov)


if __name__ == '__main__':
    main()
