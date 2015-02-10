// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "PolarGrid.h"

namespace local = turbooctospice;

local::PolarGrid::PolarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, 
lk::AbsBinningCPtr axis3) : AbsTwoPointGrid(axis1, axis2, axis3) {
    x1minSq = xmin[0]*xmin[0];
    x1maxSq = xmax[0]*xmax[0];
}

local::PolarGrid::~PolarGrid() { }
