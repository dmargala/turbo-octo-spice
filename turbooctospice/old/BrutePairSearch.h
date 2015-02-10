// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
#define TURBOOCTOSPICE_BRUTE_PAIR_SEARCH

#include <utility>

#include <iostream>
#include <iterator>

#include "types.h"

#include "boost/progress.hpp"

namespace turbooctospice {

    /// BrutePairSearch is a specific type.
    /// \tparam PixelIterable
    ///
    template <class T>
    class BrutePairSearch {
    public:
        typedef T PixelIterable;
        typedef typename PixelIterable::iterator PixelIterator;
        BrutePairSearch(PixelIterable a, PixelIterable b, bool verbose = false) : 
        _a(a), _b(b), _verbose(verbose) {
            _auto = false;
            int n(_a.size()), m(_b.size());
            _a_it = _a.begin();
            _b_it = _b.begin();
            if (_verbose) std::cout << "Number of distinct pairs : " << n*m << std::endl;
        };
        BrutePairSearch(PixelIterable a, bool verbose = false) : 
        _a(a), _verbose(verbose) {
            _auto = true;
            int n(_a.size());
            _a_it = _a.begin();
            _b_it = std::next(_a_it, 1);
            if (_verbose) std::cout << "Number of distinct pairs : " << n*(n-1)/2 << std::endl;
        };
        bool valid() const {
            if (_b_it != _a.end()) return true;
            return false;
        }
        void next() {
            ++_b_it;
            if (_b_it == _a.end()) {
                ++_a_it;
                _b_it = std::next(_a_it, 1);
            }
        }
        template <class PairType> PairType get() const {
            return PairType(*_a_it, *_b_it);
        };
    private:
        bool _verbose, _auto;
        PixelIterable _a, _b;
        PixelIterator _a_it, _b_it;
    }; // BrutePairSearch

    typedef BrutePairSearch<Pixels> BruteSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BRUTE_PAIR_SEARCH
