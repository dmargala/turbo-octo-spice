// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "BrutePairSearch.h"

namespace local = turbooctospice;

local::BrutePairSearch::BrutePairSearch() { }

local::BrutePairSearch::~BrutePairSearch() { }

void local::BrutePairSearch::findPairs(local::PairGenerator::caller_type& yield, 
std::vector<local::Pixel> const &pixels) const {
    for(int i = 0; i < pixels.size()-1; ++i) {
        for(int j = i+1; j < pixels.size(); ++j) {
            yield(local::PixelPair(pixels[i],pixels[j]));
        }
    }
}
