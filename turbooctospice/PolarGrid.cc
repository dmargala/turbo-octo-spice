// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "PolarGrid.h"

namespace local = turbooctospice;

local::PolarGrid::PolarGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, 
likely::AbsBinningCPtr axis3) : AbsTwoPointGrid(axis1, axis2, axis3) {
    x1minSq = xmin[0]*xmin[0];
    x1maxSq = xmax[0]*xmax[0];
}

local::PolarGrid::~PolarGrid() { }

bool local::PolarGrid::getSeparation(ForestPixel const &a, ForestPixel const &b,
double const &cosij, double const &thetaij, std::vector<double> &separation) const {
    double distSq = a.distance*a.distance + b.distance*b.distance - 2*a.distance*b.distance*cosij;
    if(distSq >= x1maxSq || distSq < x1minSq) return false;
    separation[0] = std::sqrt(distSq);
    separation[1] = std::fabs(a.distance-b.distance)/separation[0];
    if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
    separation[2] = 0.5*(a.wavelength+b.wavelength) - logLyA;
    if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
    return true;
}

bool local::PolarGrid::getBinIndex(ForestPixel const &a, ForestPixel const &b,
double const &cosij, double const &thetaij, int &binIndex) const {
    binIndex = 0;
    double distSq = a.distance*a.distance + b.distance*b.distance - 2*a.distance*b.distance*cosij;
    if(distSq >= x1maxSq || distSq < x1minSq) return false;
    double dist = std::sqrt(distSq);
    binIndex = int((dist - 0)/4.0);
    // double mu = std::fabs(a.distance-b.distance)/dist;
    // if(mu < xmin[1] || mu >= xmax[1]) return false;
    // binIndex = int((mu - 0)/1.0) + binIndex*50;
    // double z = 0.5*(a.wavelength+b.wavelength) - logLyA;
    // if(z < xmin[2] || z >= xmax[2]) return false;
    // binIndex = 0 + binIndex*1;
    return true;
}
