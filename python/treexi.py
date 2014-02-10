#!/usr/bin/env python

#Created 30-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
#KDTree test program for correlation function estimation


from scipy import spatial
import numpy

import argparse

def treexi(data, nbins=50., rmin=0., rmax=200.):
    print 'Starting treexi...'
    points = data[:,0:3]
    tree = spatial.KDTree(points)
    print 'Finished making tree'

    spacing = (rmax-rmin)/nbins
    dsum = numpy.zeros(nbins)
    wsum = numpy.zeros(nbins)

    nused = 0

    pairs = tree.query_pairs(rmax,eps=.01)

    npairs = len(pairs)

    print npairs
    for i,j in pairs:
        d = numpy.linalg.norm(points[i]-points[j])
        index = int((d-rmin)/spacing)
        if index < 0 or index >= nbins:
            continue
        wgt = data[i,4]*data[j,4]
        dsum[index] += wgt*data[i,3]*data[j,3]
        wsum[index] += wgt
        nused += 1;
    print "used %d of %d pairs." % (nused, npairs)
    for i in range(nbins):
        if wsum[i] > 0:
            dsum[i] /= wsum[i]

    return dsum

def main():
    """Example usage"""

    # parse command-line args
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input", type = str, default = None,
        help = "input file name")
    parser.add_argument("-r", type = float, default = 5,
        help = "max distance")
    args = parser.parse_args()

    try:
        data = numpy.loadtxt(args.input,dtype=numpy.float32)
        print 'Read %d lines from %s' % (len(data), args.input)
    except Exception,e:
        raise RuntimeError('Failed to load file %r' % args.input)

    xi = treexi(data)

    print xi

if __name__ == '__main__':
    main()