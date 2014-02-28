// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Example strategy (polcicy) design pattern using templated classes

#include "turbooctospice.h"

namespace tos = turbooctospice;

int main(int argc, char **argv) {

	// Create sample data
	std::vector<tos::Pixel> pixels;
	for(int i = 0; i < 5; ++i) {
		pixels.push_back(tos::Pixel(i));
	}

    // Example 1 
    std::cout << "Example 1: " << std::endl;
    typedef tos::XiEstimator<tos::BrutePairSearch, tos::DoNothing> SimpleXi;
 
    SimpleXi simple;
    simple.run(pixels); 
 
    // Example 2
    std::cout << "Example 2: " << std::endl;
    typedef tos::XiEstimator<tos::BrutePairSearch, tos::PrintPair> WeightedXi;
 
    WeightedXi weighted;
    weighted.run(pixels);
}