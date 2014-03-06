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
    template <
        class PixelType, class PixelIterable,
        template <class> class PairSearchPolicy, 
        template <class> class BinPolicy
    > 
    class XiEstimator {
    public:
        typedef std::pair<PixelType,PixelType> PixelPair;
        typedef boost::coroutines::coroutine<PixelPair()> PairGenerator;
        /// Creates a new algorithm templated on the specified PairSearchPolicy and BinPolicy.
        /// Set the verbose option to true to print pair statistics to console.
        /// @param psp a PairSearchPolicy pointer
        /// @param bp a BinPolicy pointer
        /// @param verbose a boolean argument.
        ///
        XiEstimator(PairSearchPolicy<PixelIterable> *psp, BinPolicy<PixelType> *bp, bool verbose = false) : 
        _psp(psp), _bp(bp), _verbose(verbose) {};
        void run(PixelIterable const &a, PixelIterable const &b, std::vector<double> &xi) const {
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
            for(int index = 0; index < nbins; ++index) {
                if(wsum[index] > 0) dsum[index] /= wsum[index];
            }
            dsum.swap(xi);
        }
        PairGenerator getGenerator(PixelIterable const &a, PixelIterable const &b) const {
            if(&a == &b) {
                int n(a.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*(n-1)/2 << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy<PixelIterable>::template findPairs<PairGenerator>, _psp, _1, a));
            }
            else {
                int n(a.size()), m(b.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*m << std::endl;
                return PairGenerator(boost::bind(&PairSearchPolicy<PixelIterable>::template findPairs<PairGenerator>, _psp, _1, a, b));
            }
        }
    private:
        boost::shared_ptr<const PairSearchPolicy<PixelIterable> > _psp;     ///< PairSearchPolicy pointer
        boost::shared_ptr<const BinPolicy<PixelType> > _bp;             ///< BinPolicy pointer
        bool _verbose;                                      ///< verbose

    }; // XiEstimator
    
    typedef XiEstimator<Pixel, Pixels, BrutePairSearch, IgnorePair> BruteIgnoreXi;
    typedef XiEstimator<Pixel, Pixels, BrutePairSearch, BinXYZPair> BruteBinXi;
    typedef XiEstimator<Pixel, Pixels, BucketPairSearch, IgnorePair> BucketIgnoreXi;
    typedef XiEstimator<Pixel, Pixels, BucketPairSearch, BinXYZPair> BucketBinXi; 

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR