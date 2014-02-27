#include "XiEstimator.h"

#include <iostream>
#include <string>

namespace local = turbooctospice;

float local::BinPolicyDummy::binPair(float a, float b) {
    return a*b;
}

float local::BinPolicyWeighted::binPair(float a, float b) {
    float weight = .5;
    return a*b*weight*weight;
}
