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
    class XiEstimator {
    public:
        XiEstimator(PairSearchPolicy *psp, BinPolicy *bp) : _psp(psp), _bp(bp) {};
        void run(Pixels const &a, Pixels const &b, std::vector<double> &xi) const {
            // create internal accumulation vectors
            int nbins = _bp->getNBinsTotal();
            std::vector<double> dsum(nbins,0.), wsum(nbins,0.);

            // Loop over pixel pairs
            PairGenerator pairs = getGenerator(a, b);
            long npair(0), nused(0);
            while(pairs){ // Check completion status
                npair++;
                // Extract yielded result
                _bp->binPair(pairs.get().first, pairs.get().second, dsum, wsum, nused);
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
                int n(a.size());
                std::cout << "Number of distinct pairs : " << n*(n-1)/2 << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy::findPairs, _psp, _1, a));
            }
            else {
                int n(a.size()), m(b.size());
                std::cout << "Number of distinct pairs : " << n*m << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy::findPairs, _psp, _1, a, b));
            }
        }
    private:
        boost::shared_ptr<const PairSearchPolicy> _psp;
        boost::shared_ptr<const BinPolicy> _bp;

    }; // XiEstimator

    typedef XiEstimator<BruteSearch, Ignore> BruteIgnoreXi;
    typedef XiEstimator<BruteSearch, Bin> BruteBinXi;
    typedef XiEstimator<BucketSearch, Ignore> BucketIgnoreXi;
    typedef XiEstimator<BucketSearch, Bin> BucketBinXi; 

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR