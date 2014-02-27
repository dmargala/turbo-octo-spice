#include "XiEstimator.h"

#include <iostream>
#include <string>

namespace local = turbooctospice;

// // Behaviour method
// template <typename PairSearchPolicy, typename BinPolicy> void local::XiEstimator::run() const {
//     // Two policy methods
//     print(message());
// };

// template<typename MessageType> void local::PairSearchPolicyBrute::print(MessageType const &message) const {
//     std::cout << message << std::endl;
// }

std::string local::BinPolicyDummy::message() const {
    return "Hello, World!";
}

std::string local::BinPolicyWeighted::message() const {
    return "Hello, Weighted World!";
}
