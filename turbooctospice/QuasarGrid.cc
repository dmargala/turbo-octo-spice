// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "QuasarGrid.h"

namespace local = turbooctospice;

local::QuasarGrid::QuasarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, 
lk::AbsBinningCPtr axis3) : AbsTwoPointGrid(axis1, axis2, axis3) { }

local::QuasarGrid::~QuasarGrid() { }
