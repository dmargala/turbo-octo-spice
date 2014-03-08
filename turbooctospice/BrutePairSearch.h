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
		BrutePairSearch(PixelIterable a, PixelIterable b, bool verbose = false) : 
		_a(a), _b(b), _verbose(verbose) {
			_auto = false;
		};
		BrutePairSearch(PixelIterable a, bool verbose = false) : 
		_a(a), _verbose(verbose) {
			_auto = true;
		};
        template <class PairGenerator, class PairType> void findPairs(typename PairGenerator::caller_type& yield) const {
			if (_auto) {
				if (_verbose) std::cout << "Entering auto-correlation generator ..." << std::endl;
				int n(_a.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*(n-1)/2 << std::endl;
	    		for(auto i = _a.begin(); i != boost::prior(_a.end()); ++i) {
			        for(auto j = boost::next(i); j != _a.end(); ++j) {
			            yield(PairType(*i,*j));
			        }
			    }
			    if (_verbose) std::cout << "Exiting auto-correlation generator ..." << std::endl;
			}
			else {
	        	if (_verbose) std::cout << "Entering cross-correlation generator ..." << std::endl;
	        	int n(_a.size()), m(_b.size());
                if (_verbose) std::cout << "Number of distinct pairs : " << n*m << std::endl;
	    		for(auto i = _a.begin(); i != _a.end(); ++i) {
			        for(auto j = _b.begin(); j != _b.end(); ++j) {
			            yield(PairType(*i,*j));
			        }
			    }
			    if (_verbose) std::cout << "Exiting cross-correlation generator ..." << std::endl;
		   	}
        }
	private:
		bool _verbose, _auto;
		PixelIterable _a, _b;
	}; // BrutePairSearch

	typedef BrutePairSearch<Pixels> BruteSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
