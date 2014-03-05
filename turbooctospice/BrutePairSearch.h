// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#pragma clang diagnostic ignored "-Wc++11-extensions"

#ifndef TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
#define TURBOOCTOSPICE_BRUTE_PAIR_SEARCH

#include <utility>

#include <iostream>

#include "types.h"

namespace turbooctospice {

	template <typename PixelIterable>
	class BrutePairSearch : private PixelIterable {
	public:
        void findPairs(PairGenerator::caller_type& yield, PixelIterable const &a, PixelIterable const &b) const {
        	std::cout << "Entering cross-correlation generator ..." << std::endl;
    		for(auto i = a.begin(); i != a.end(); ++i) {
		        for(auto j = b.begin(); j != b.end(); ++j) {
		            yield(std::make_pair(*i,*j));
		        }
		    }
		    std::cout << "Exiting cross-correlation generator ..." << std::endl;
        }
        void findPairs(PairGenerator::caller_type& yield, PixelIterable const &a) const {
        	std::cout << "Entering auto-correlation generator ..." << std::endl;
    		for(auto i = a.begin(); i != boost::prior(a.end()); ++i) {
		        for(auto j = boost::next(i); j != a.end(); ++j) {
		            yield(std::make_pair(*i,*j));
		        }
		    }
		    std::cout << "Exiting auto-correlation generator ..." << std::endl;
        }
	private:
	}; // BrutePairSearch

	typedef BrutePairSearch<Pixels> BruteSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
