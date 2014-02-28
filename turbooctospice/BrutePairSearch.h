// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
#define TURBOOCTOSPICE_BRUTE_PAIR_SEARCH

#include <vector>

#include "types.h"

namespace turbooctospice {
	class BrutePairSearch {
	public:
		BrutePairSearch();
		virtual ~BrutePairSearch();
        void findPairs(PairGenerator::caller_type& yield, std::vector<Pixel> const &pixels) const;
	private:
	}; // BrutePairSearch
} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
