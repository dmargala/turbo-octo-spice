// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TOS_XI_ESTIMATOR
#define TOS_XI_ESTIMATOR

#include <iostream>
#include <string>
#include <vector>

#include "boost/bind.hpp"

#include "types.h"

namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : private PairSearchPolicy, private BinPolicy {
        // using PairSearchPolicy::findPairs;
        // using BinPolicy::binPair;
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
    };

    class PairSearchPolicyBrute {
    public:
        void findPairs(PairGenerator::caller_type& yield, std::vector<Pixel> const &pixels) const;
    };


    class BinPolicyDummy {
    protected:
        void binPair(PixelPair const &pair) const;
    };

    class BinPolicyWeighted {
    protected:
        void binPair(PixelPair const &pair) const; 
    };

}

#endif