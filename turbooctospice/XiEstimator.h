// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include "boost/bind.hpp"

#include "types.h"

namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : private PairSearchPolicy, private BinPolicy {
    public:
        void run(Pixels const &a, Pixels const &b) {
            PairGenerator pairs(boost::bind(&PairSearchPolicy::findPairs, 
                dynamic_cast<PairSearchPolicy*>(this), _1, a, b));
            // Loop over yielded pairs
            while(pairs){ // Check completion status
                // Extract yielded result
                BinPolicy::binPair(pairs.get());
                // Fire up the generator?
                pairs();
            }
        }
        void run(Pixels const &a) {
            run(a, a);
        }
    private:
    }; // XiEstimator

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR