#include "tos.h"

namespace tos = turbooctospice;

int main(int argc, char **argv) {
    /* Example 1 */
    typedef tos::XiEstimator<tos::PairSearchPolicyBrute, tos::BinPolicyDummy> SimpleXi;
 
    SimpleXi simple;
    simple.run(); // prints "Hello, World!"
 
    /* Example 2 */
    typedef tos::XiEstimator<tos::PairSearchPolicyBrute, tos::BinPolicyWeighted> WeightedXi;
 
    WeightedXi weighted;
    weighted.run(); // prints "Hello, Weighted Worled!"
}