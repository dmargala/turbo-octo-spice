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
    /// @tparam PairSearchPolicy policy class must expose a findPairs member function
    /// @tparam BinPolicy policy class must expose an inner type called PairType and a binPair member function
    ///
    template <class PairSearchPolicy, class BinPolicy> 
    class XiEstimator {
        typedef typename BinPolicy::PairType PairType;
        typedef typename PairType::PairGenerator PairGenerator;
    public:
        /// Creates a new algorithm templated on the specified PairSearchPolicy and BinPolicy.
        /// Set the verbose option to true to print pair statistics to console.
        /// @param psp a PairSearchPolicy pointer
        /// @param bp a BinPolicy pointer
        /// @param verbose a boolean argument.
        ///
        XiEstimator(PairSearchPolicy *psp, BinPolicy *bp, bool verbose = false) : _psp(psp), _bp(bp), _verbose(verbose) {};

        /// Perform xi estimate
        /// @param a iterable container of pixel data
        /// @param b iterable container of pixel data
        /// @param xi output xi estimate
        /// @param normalize divide bin totals by weights (turn off to count pairs as a function of distance)
        ///
        std::vector<double> run(bool normalize = true) const {
            // create internal accumulation vectors
            int nbins = _bp->getNBins();
            std::vector<double> dsum(nbins,0.), wsum(nbins,0.);

            // Loop over pixel pairs
            // Binds the PairSearchPolicy findPairs method to a generator of PixelPairs
            PairGenerator pairs(boost::bind(&PairSearchPolicy::template findPairs<PairGenerator, PairType>, _psp, _1));

            long npair(0), nused(0);
            while(pairs){ // Check completion status
            //for(auto &pair : pairs) {
                npair++;
                // Extract yielded result
                _bp->binPair(pairs.get(), dsum, wsum, nused);
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
            return dsum;
        }
    private:
        boost::shared_ptr<const PairSearchPolicy> _psp;    
        boost::shared_ptr<const BinPolicy> _bp;     
        bool _verbose;             

    }; // XiEstimator
    
    typedef XiEstimator<BruteSearch, Ignore> BruteIgnoreXi;
    typedef XiEstimator<BruteSearch, BinXYZ> BruteBinXi;
    typedef XiEstimator<BucketSearch, Ignore> BucketIgnoreXi;
    typedef XiEstimator<BucketSearch, BinXYZ> BucketBinXi;

#ifdef HAVE_LIBHEAL
    typedef XiEstimator<HealSearch, IgnoreAng > HealIgnoreXi;
    typedef XiEstimator<HealSearch, BinAng> HealBinXi;
#endif

} // turbooctospice 

#endif // TURBOOCTOSPICE_XI_ESTIMATOR