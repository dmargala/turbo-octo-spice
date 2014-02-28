// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <vector>

#include "boost/bind.hpp"

#include "types.h"

namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : public PairSearchPolicy, public BinPolicy {
    public:
        void run(std::vector<Pixel> const &pixels) {
            PairGenerator pairs(boost::bind(&PairSearchPolicy::findPairs, 
                dynamic_cast<PairSearchPolicy*>(this), _1, pixels));
            // Loop over yielded pairs
            while(pairs){ // Check completion status
                // Extract yielded result
                BinPolicy::binPair(pairs.get());
                // Fire up the generator?
                pairs();
            }
        }
    private:
    }; // XiEstimator

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR