// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "AbsTwoPointGrid.h"

namespace local = turbooctospice;


local::AbsTwoPointGrid::AbsTwoPointGrid(lk::AbsBinningCPtr axis1, 
lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3) : _grid(axis1, axis2, axis3) {
    for(int axis = 0; axis < 3; ++axis) {
        auto bins = _grid.getAxisBinning(axis);
        xmin.push_back(bins->getBinLowEdge(0));
        xmax.push_back(bins->getBinHighEdge(bins->getNBins()-1));
    }
}

local::AbsTwoPointGrid::~AbsTwoPointGrid() { }

