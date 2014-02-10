#!/usr/bin/env python

#Created 31-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

import argparse
import numpy

lymanAlpha = 1216.

def main():
    """Add fake line of sight data to qso position list"""

    # parse command-line args
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i","--input", type = str, default = None,
        help = "input file name")
    parser.add_argument("-o","--output", type = str, default = None,
        help = "output file name")
    parser.add_argument("--forest-high", type = float, default = 1200.,
        help = "Lyman-alpha forest high cutoff wavelength")
    parser.add_argument("--forest-low", type = float, default = 1050.,
        help = "Lyman-alpha forest low cutoff wavelength")
    parser.add_argument("--spec-low", type = float, default = 3650.,
        help = "Spectrograph wavelength lower-limit")
    parser.add_argument("--combine", type = int, default = 3,
        help = "Number of observed wavelengths bins to combine.")
    args = parser.parse_args()

    try:
        data = numpy.loadtxt(args.input,dtype=numpy.float32)
        print 'Read %d lines from %s' % (len(data), args.input)
    except Exception,e:
        raise RuntimeError('Failed to load file %r' % args.input)

    m = (10**(.0001))**args.combine

    out = open(args.output, 'w')

    qsocounter = 0
    counter = 0
    pixels = []
    for quasar in data:
        qsocounter += 1
        if (qsocounter % 10000) == 0:
            print 'Processing quasar number %d, written %d pixels so far... ' % (qsocounter, counter)
        zquasar = quasar[2]
        if zquasar < 2:
            continue
        zmin = max(args.spec_low/lymanAlpha-1, args.forest_low/lymanAlpha*(1+zquasar)-1)
        zmax = args.forest_high/lymanAlpha*(1+zquasar)-1
        z = zmin
        while z < zmax:
            z = (1+z)*m-1
            out.write('%f %f %f\n' % (quasar[0], quasar[1], z))
            counter += 1
    out.close()
    print 'Wrote %d pixels to %s' % (counter,args.output)

if __name__ == '__main__':
    main()