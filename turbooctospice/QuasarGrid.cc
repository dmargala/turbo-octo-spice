// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "QuasarGrid.h"

namespace local = turbooctospice;

local::QuasarGrid::QuasarGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
likely::AbsBinningCPtr axis3) : AbsTwoPointGrid(axis1, axis2, axis3) { }

local::QuasarGrid::~QuasarGrid() { }

bool local::QuasarGrid::getSeparation(ForestPixel const &a, ForestPixel const &b,
double const &cosij, double const &thetaij, std::vector<double> &separation) const {
    separation[0] = std::fabs(a.loglam-b.loglam);
    if(separation[0] >= xmax[0] || separation[0] < xmin[0]) return false;
    separation[1] = thetaij*rad2arcmin;
    // we don't need to check thetaij, we've already done this for the line of sights
    // if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
    separation[2] = 0.5*(a.loglam+b.loglam) - logLyA;
    if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
    return true;
}

bool local::QuasarGrid::getBinIndex(ForestPixel const &a, ForestPixel const &b,
double const &cosij, double const &thetaij, int &binIndex) const {
    return false;
}
