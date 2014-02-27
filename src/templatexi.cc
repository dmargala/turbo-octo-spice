#include "tos.h"

namespace tos = turbooctospice;

int main(int argc, char **argv) {

	std::vector<float> pixels;
	for(int i = 0; i < 5; ++i) {
		pixels.push_back(i);
	}

    /* Example 1 */
    std::cout << "Example 1: " << std::endl;
    typedef tos::XiEstimator<tos::PairSearchPolicyBrute, tos::BinPolicyDummy> SimpleXi;
 
    SimpleXi simple;
    simple.run(pixels); 
 
    /* Example 2 */
    std::cout << "Example 2: " << std::endl;
    typedef tos::XiEstimator<tos::PairSearchPolicyBrute, tos::BinPolicyWeighted> WeightedXi;
 
    WeightedXi weighted;
    weighted.run(pixels);
}