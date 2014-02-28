// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "XiEstimator.h"

#include <iostream>
#include <string>

namespace local = turbooctospice;

// template<typename BinMethod> void local::PairSearchPolicyBrute::findPairs(std::vector<float> &pixels, BinMethod binPair) const {
//     for(int i = 0; i < pixels.size()-1; ++i) {
//         float a = pixels[i];
//         for(int j = i+1; j < pixels.size(); ++j) {
//             float b = pixels[j];
//             std::cout << a << "," << b << " -> " << binPair(a, b) << std::endl;
//         }
//     }
// }

void local::BinPolicyDummy::binPair(local::PixelPair const &pair) const {
    std::cout << pair.first*pair.second << std::endl;
}

void local::BinPolicyWeighted::binPair(local::PixelPair const &pair) const {
    float weight = .5;
    std::cout << pair.first*pair.second*weight*weight << std::endl;
}
