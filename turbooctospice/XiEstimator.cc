// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "XiEstimator.h"

#include <iostream>
#include <string>

namespace local = turbooctospice;

void local::PairSearchPolicyBrute::findPairs(local::PairGenerator::caller_type& yield, 
std::vector<local::Pixel> const &pixels) const {
    for(int i = 0; i < pixels.size()-1; ++i) {
        for(int j = i+1; j < pixels.size(); ++j) {
            yield(local::PixelPair(pixels[i],pixels[j]));
        }
    }
}

void local::BinPolicyDummy::binPair(local::PixelPair const &pair) const {
    std::cout << pair.first*pair.second << std::endl;
}

void local::BinPolicyWeighted::binPair(local::PixelPair const &pair) const {
    float weight = .5;
    std::cout << pair.first*pair.second*weight*weight << std::endl;
}
