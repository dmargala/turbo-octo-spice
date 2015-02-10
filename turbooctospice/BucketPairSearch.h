// Created 04-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
#define TURBOOCTOSPICE_BUCKET_PAIR_SEARCH

#include "likely/likely.h"

#include "boost/foreach.hpp"
#include "boost/progress.hpp"

#include <iostream>

namespace lk = likely;

namespace turbooctospice {

    template <class T>
    class BucketPairSearch {
    public:
        typedef T PixelIterable;
        typedef typename PixelIterable::iterator PixelIterator;
        typedef std::map<int, std::vector<int> > BucketToPoints;
        typedef std::map<int, std::vector<int> > BucketToBuckets;
        typedef typename std::vector<int>::iterator PointIterator;
        typedef typename std::vector<int>::iterator BucketIterator;
        typedef typename BucketToPoints::iterator BucketMapIterator;
        BucketPairSearch(PixelIterable a, PixelIterable b, lk::BinnedGrid const &bucketgrid, bool verbose = false) : 
        _a(a), _b(b), _bucketgrid(bucketgrid), _verbose(verbose) {
            _auto = false;
            initialize();
        };
        BucketPairSearch(PixelIterable a, lk::BinnedGrid const &bucketgrid, bool verbose = false) : 
        _a(a), _bucketgrid(bucketgrid), _verbose(verbose) { 
            _auto = true;
            initialize();
        };
        ~BucketPairSearch() {};
        void initialize() {
            std::vector<double> position(3);
            std::vector<int> binNeighbors;
            for(int i = 0; i < _a.size(); ++i) {
                position[0] = _a.at(i).x;
                position[1] = _a.at(i).y;
                position[2] = _a.at(i).z;
                int bucketIndex = _bucketgrid.getIndex(position);
                if(_bucketPointsMapA.count(bucketIndex) > 0) {
                    _bucketPointsMapA[bucketIndex].push_back(i);
                }
                else {
                    _bucketPointsMapA[bucketIndex] = std::vector<int>(1,i);
                    _bucketgrid.getBinNeighbors(bucketIndex, binNeighbors);
                    _bucketNeighborsMap[bucketIndex] = binNeighbors;
                }
            }

            if (_auto) {
                _bucketPointsMapB = _bucketPointsMapA;
            }
            else {
                for(int i = 0; i < _b.size(); ++i) {
                    position[0] = _b.at(i).x;
                    position[1] = _b.at(i).y;
                    position[2] = _b.at(i).z;
                    int bucketIndex = _bucketgrid.getIndex(position);
                    if(_bucketPointsMapB.count(bucketIndex) > 0) {
                        _bucketPointsMapB[bucketIndex].push_back(i);
                    }
                    else {
                        _bucketPointsMapB[bucketIndex] = std::vector<int>(1,i);
                        _bucketgrid.getBinNeighbors(bucketIndex, binNeighbors);
                        _bucketNeighborsMap[bucketIndex] = binNeighbors;
                    }
                }

            }   

            int nbuckets = _bucketNeighborsMap.size();
            std::cout << "We have " << nbuckets << " buckets" << std::endl;

            _aBucketIt = _bucketPointsMapA.begin();
            _aPointsIt = (_aBucketIt->second).begin();
            _neighboringBucketsIt = (_bucketNeighborsMap[_aBucketIt->first]).begin();
            _bPointsIt = (_bucketPointsMapB[*_neighboringBucketsIt]).begin();
        };
        bool valid() const {
            if (_aBucketIt != _bucketPointsMapA.end()) return true;
            return false;
        }
        void next() {
            ++_bPointsIt; //int j
            if (_bPointsIt == (_bucketPointsMapB[*_neighboringBucketsIt]).end()) {
                ++_neighboringBucketsIt; //int bucketindex
                if (_neighboringBucketsIt == (_bucketNeighborsMap[_aBucketIt->first]).end()) {
                    ++_aPointsIt; //int i 
                    if (_aPointsIt == (_aBucketIt->second).end()) {
                        ++_aBucketIt; //auto bucket
                        _aPointsIt = (_aBucketIt->second).begin();
                    }
                    _neighboringBucketsIt = (_bucketNeighborsMap[_aBucketIt->first]).begin();
                }
                if(_auto && (*_bPointsIt <= *_aPointsIt)) next();
                _bPointsIt = (_bucketPointsMapB[*_neighboringBucketsIt]).begin();
            }
        }
        template <class PairType> PairType get() const {
            return PairType(_a.at(*_aPointsIt), _b.at(*_bPointsIt));
        };
    private:
        bool _verbose, _auto;
        PixelIterable _a, _b;
        PointIterator _aPointsIt, _bPointsIt;
        BucketIterator _neighboringBucketsIt;
        BucketMapIterator _aBucketIt;
        // The key is a global bucketgrid index and the value is a 
        // list of indices that represent points inside that bucket
        BucketToPoints _bucketPointsMapA, _bucketPointsMapB;
        // The key is a global bucketgrid index and the value is
        // a list of neighboring buckets
        BucketToBuckets _bucketNeighborsMap;
        lk::BinnedGrid _bucketgrid;

    }; // BucketPairSearch

    typedef BucketPairSearch<Pixels> BucketSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
