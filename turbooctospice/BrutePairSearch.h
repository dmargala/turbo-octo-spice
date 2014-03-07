// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#pragma clang diagnostic ignored "-Wc++11-extensions"

#ifndef TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
#define TURBOOCTOSPICE_BRUTE_PAIR_SEARCH

#include <utility>

#include <iostream>

#include "types.h"

namespace turbooctospice {

	/// BrutePairSearch is a specific type.
	/// \tparam PixelIterable
	///
	template <class T>
	class BrutePairSearch {
	public:
		typedef T PixelIterable;
		BrutePairSearch(bool verbose = false) : _verbose(verbose) {};
        template <class PairGenerator> void findPairs(typename PairGenerator::caller_type& yield, PixelIterable const &a, PixelIterable const &b) const {
        	if (_verbose) std::cout << "Entering cross-correlation generator ..." << std::endl;
    		for(auto i = a.begin(); i != a.end(); ++i) {
		        for(auto j = b.begin(); j != b.end(); ++j) {
		            yield(std::make_pair(*i,*j));
		        }
		    }
		    if (_verbose) std::cout << "Exiting cross-correlation generator ..." << std::endl;
        }
        template <class PairGenerator> void findPairs(typename PairGenerator::caller_type& yield, PixelIterable const &a) const {
        	if (_verbose) std::cout << "Entering auto-correlation generator ..." << std::endl;
    		for(auto i = a.begin(); i != boost::prior(a.end()); ++i) {
		        for(auto j = boost::next(i); j != a.end(); ++j) {
		            yield(std::make_pair(*i,*j));
		        }
		    }
		    if (_verbose) std::cout << "Exiting auto-correlation generator ..." << std::endl;
        }
	private:
		bool _verbose;
	}; // BrutePairSearch

	typedef BrutePairSearch<Pixels> BruteSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
