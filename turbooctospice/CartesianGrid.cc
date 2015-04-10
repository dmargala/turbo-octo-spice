// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "CartesianGrid.h"

namespace local = turbooctospice;

local::CartesianGrid::CartesianGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, 
likely::AbsBinningCPtr axis3) : AbsTwoPointGrid(axis1, axis2, axis3) {
    x1minSq = xmin[0]*xmin[0];
    x1maxSq = xmax[0]*xmax[0];
};


local::CartesianGrid::~CartesianGrid() { }

bool local::CartesianGrid::getSeparation(ForestPixel const &a, ForestPixel const &b, 
double const &cosij, double const &thetaij, std::vector<double> &separation) const {
    double distSq = a.distance*a.distance + b.distance*b.distance - 2*a.distance*b.distance*cosij;
    if(distSq >= x1maxSq || distSq < x1minSq) return false;
    separation[0] = std::fabs(a.distance-b.distance);
    separation[1] = std::sqrt(distSq - separation[0]*separation[0]);
    // separation[1] = thetaij*cosmology->getTransverseComovingScale(separation[2]);
    if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
    separation[2] = 0.5*(a.loglam+b.loglam) - logLyA;
    if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
    return true;
}


bool local::CartesianGrid::getBinIndex(ForestPixel const &a, ForestPixel const &b,
double const &cosij, double const &thetaij, int &binIndex) const {
    return false;
}