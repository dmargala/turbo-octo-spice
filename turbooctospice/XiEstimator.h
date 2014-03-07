// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <iostream> 

#include "boost/bind.hpp"

#include "types.h"

#include "likely/likely.h"

namespace lk = likely;

namespace turbooctospice {
    /// Represents a templated algorithm for estimating a correlation function from a set of pixels.
    /// \tparam PairSearchPolicy psp
    /// \tparam BinPolicy bp
    ///
    template <typename PairSearchPolicy, typename BinPolicy> 
    class XiEstimator {
        typedef typename PairSearchPolicy::PixelIterable PixelIterable;
        typedef typename BinPolicy::PixelType PixelType;
        typedef std::pair<PixelType,PixelType> PixelPair;
        typedef boost::coroutines::coroutine<PixelPair()> PairGenerator;
    public:
        /// Creates a new algorithm templated on the specified PairSearchPolicy and BinPolicy.
        /// Set the verbose option to true to print pair statistics to console.
        /// @param psp a PairSearchPolicy pointer
        /// @param bp a BinPolicy pointer
        /// @param verbose a boolean argument.
        ///
        XiEstimator(PairSearchPolicy *psp, BinPolicy *bp, bool verbose = false) : 
        _psp(psp), _bp(bp), _verbose(verbose) {};
        void run(PixelIterable const &a, PixelIterable const &b, std::vector<double> &xi, bool normalize = true) const {
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

            if (_verbose) std::cout << "used " << nused << " of " << npair << " pairs." << std::endl;
            // Compute xi and swap result with xi reference argument
            if (normalize) {
                for(int index = 0; index < nbins; ++index) {
                    if(wsum[index] > 0) dsum[index] /= wsum[index];
                }
            }
            dsum.swap(xi);
        }
        PairGenerator getGenerator(PixelIterable const &a, PixelIterable const &b) const {
            if(&a == &b) {
                int n(a.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*(n-1)/2 << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy::template findPairs<PairGenerator>, _psp, _1, a));
            }
            else {
                int n(a.size()), m(b.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*m << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy::template findPairs<PairGenerator>, _psp, _1, a, b));
            }
        }
    private:
        boost::shared_ptr<const PairSearchPolicy> _psp;     ///< PairSearchPolicy pointer
        boost::shared_ptr<const BinPolicy> _bp;             ///< BinPolicy pointer
        bool _verbose;                                      ///< verbose

    }; // XiEstimator
    
    typedef XiEstimator<BruteSearch, Ignore> BruteIgnoreXi;
    typedef XiEstimator<BruteSearch, Bin> BruteBinXi;
    typedef XiEstimator<BucketSearch, Ignore> BucketIgnoreXi;
    typedef XiEstimator<BucketSearch, Bin> BucketBinXi; 

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR