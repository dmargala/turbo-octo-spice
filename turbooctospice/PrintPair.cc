// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "PrintPair.h"

#include <iostream>

namespace local = turbooctospice;

local::PrintPair::PrintPair() { }

local::PrintPair::~PrintPair() { }

void local::PrintPair::binPair(PixelPair const &pair) const {
	std::cout << pair.first << "," << pair.second << std::endl;
};
