// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_CONSTANTS
#define TURBOOCTOSPICE_CONSTANTS

#include <cmath>

namespace turbooctospice {

    const double lyA = 1216.0;
    const double logLyA = std::log10(lyA);
    const double pi = std::atan(1.0)*4;
    const double deg2rad = pi/180.0;
    const double rad2arcmin = 60.0/deg2rad;

}

#endif