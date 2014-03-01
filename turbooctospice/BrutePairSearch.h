// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
#define TURBOOCTOSPICE_BRUTE_PAIR_SEARCH

#include "types.h"

namespace turbooctospice {
	template <typename PixelIterable>
	class BrutePairSearch : private PixelIterable {
	public:
        void findPairs(PairGenerator::caller_type& yield, PixelIterable const &a, PixelIterable const &b) const {
    		for(typename PixelIterable::const_iterator it = a.begin(); it != a.end();++it) {
		        for(typename PixelIterable::const_iterator jt = b.begin(); jt != b.end(); ++jt) {
		            yield(PixelPair(*it,*jt));
		        }
		    }
        }
	private:
	}; // BrutePairSearch

} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
