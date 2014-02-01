#!/usr/bin/env python

#Created 31-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

import argparse
import numpy


def main():
    """Add fake line of sight data to qso position list"""

    # parse command-line args
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input", type = str, default = None,
        help = "input file name")
    parser.add_argument("-o","--output", type = str, default = None,
        help = "output file name")
    args = parser.parse_args()

    try:
        data = numpy.loadtxt(args.input,dtype=numpy.float32)
        print 'Read %d lines from %s' % (len(data), args.input)
    except Exception,e:
        raise RuntimeError('Failed to load file %r' % args.input)

    lambdaObsMin = 3650.
    lambdaAlphaCutoff = 1200.
    lambdaBetaCutoff = 1050.

    lambdaAlpha = 1216.

    pixels = []
    m = (10**(.0001))**3

    out = open(args.output, 'w')

    qsocounter = 0
    counter = 0

    for quasar in data:
        qsocounter += 1
        if (qsocounter % 10000) == 0:
            print 'Processing quasar number %d, written %d pixels so far... ' % (qsocounter, counter)
        zquasar = quasar[2]
        if zquasar < 2:
            continue
        zmin = max(lambdaObsMin/lambdaAlpha-1, lambdaBetaCutoff/lambdaAlpha*(1+zquasar)-1)
        zmax = lambdaAlphaCutoff/lambdaAlpha*(1+zquasar)-1
        z = zmin
        while z < zmax:
            z = (1+z)*m-1
            out.write('%f %f %f\n' % (quasar[0], quasar[1], z))
            counter += 1


    out.close()
    print 'Wrote %d pixels to %s' % (counter,args.output)

if __name__ == '__main__':
    main()