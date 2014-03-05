// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <iostream> 

#include "boost/bind.hpp"

#include "types.h"

#include "likely/likely.h"

namespace lk = likely;


namespace turbooctospice {
 
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator : private PairSearchPolicy, private BinPolicy {
    public:
        void run(Pixels const &a, Pixels const &b, lk::BinnedGrid const &grid, bool rmu, 
        double x1min, double x1max, double x2min, double x2max, std::vector<double> &xi) {
            // create internal accumulation vectors
            int nbins = grid.getNBinsTotal();
            std::vector<double> dsum(nbins,0.), wsum(nbins,0.);

            // Loop over pixel pairs
            PairGenerator pairs = getGenerator(a, b);
            long npair(0), nused(0);
            while(pairs){ // Check completion status
                npair++;
                // Extract yielded result
                BinPolicy::binPair(pairs.get().first, pairs.get().second, grid, rmu, 
                    x1min, x1max, x2min, x2max, dsum, wsum, nused);
                // Fire up the generator?
                pairs();
            }

            std::cout << "used " << nused << " of " << npair << " pairs." << std::endl;
            // Compute xi and swap result with xi reference argument
            for(int index = 0; index < nbins; ++index) {
                if(wsum[index] > 0) dsum[index] /= wsum[index];
            }
            dsum.swap(xi);
        }
        PairGenerator getGenerator(Pixels const &a, Pixels const &b) const {
            if(&a == &b) {
                return PairGenerator(boost::bind(&PairSearchPolicy::findPairs, this, _1, a));
            }
            else {
                return PairGenerator(boost::bind(&PairSearchPolicy::findPairs, this, _1, a, b));
            }
        }
    private:
    }; // XiEstimator

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR