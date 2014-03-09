// Created 04-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
#define TURBOOCTOSPICE_BUCKET_PAIR_SEARCH

#include "likely/likely.h"

#include "boost/foreach.hpp"

#include <iostream>


namespace lk = likely;

namespace turbooctospice {

    template <class T>
    class BucketPairSearch {
    public:
        typedef T PixelIterable;
        typedef std::map<int, std::vector<int> > BucketToPoints;
        typedef std::map<int, std::vector<int> > BucketToBuckets;
        BucketPairSearch(PixelIterable a, PixelIterable b, lk::BinnedGrid const &bucketgrid, bool verbose = false) : 
        _a(a), _b(b), _bucketgrid(bucketgrid), _verbose(verbose) {
            _auto = false;
        };
        BucketPairSearch(PixelIterable a, lk::BinnedGrid const &bucketgrid, bool verbose = false) : 
        _a(a), _bucketgrid(bucketgrid), _verbose(verbose) { 
            _auto = true;
        };
        ~BucketPairSearch() {};
        template<class PairGenerator, class PairType> void findPairs(typename PairGenerator::caller_type& yield) const {
            if (_auto) {
                if (_verbose) std::cout << "Entering auto-correlation generator ..." << std::endl;
                // The key is a global bucketgrid index and the value is a 
                // list of indices that represent points inside that bucket
                BucketToPoints bucketPointsMap;
                // The key is a global bucketgrid index and the value is
                // a list of neighboring buckets
                BucketToBuckets bucketNeighborsMap;
                // First pass through the data, fill buckets with point indices.
                // Also create a lookup tables for points->buckets and buckets->neighbors
                std::vector<double> position(3);
                std::vector<int> binNeighbors;
                for(int i = 0; i < _a.size(); ++i) {
                    position[0] = _a.at(i).x;
                    position[1] = _a.at(i).y;
                    position[2] = _a.at(i).z;
                    int bucketIndex = _bucketgrid.getIndex(position);
                    if(bucketPointsMap.count(bucketIndex) > 0) {
                        bucketPointsMap[bucketIndex].push_back(i);
                    }
                    else {
                        bucketPointsMap[bucketIndex] = std::vector<int>(1,i);
                        _bucketgrid.getBinNeighbors(bucketIndex, binNeighbors);
                        bucketNeighborsMap[bucketIndex] = binNeighbors;
                    }
                }
                int nbuckets = bucketNeighborsMap.size();
                std::cout << "We have " << nbuckets << " buckets" << std::endl;
                // Loop over all buckets
                for(auto &bucket : bucketPointsMap){
                    // Loop over all points in each bucket
                    for(int i : bucket.second) {
                        // Compare this points to all points in neighboring buckets
                        for(int bucketindex : bucketNeighborsMap[bucket.first]) {
                            // Loop over all points in neighboring bucket
                            for(int j : bucketPointsMap[bucketindex]) {
                                // Only count pairs once
                                if(j <= i) continue;
                                yield(PairType(_a.at(i),_a.at(j)));
                            }
                        }
                    }
                }
                if (_verbose) std::cout << "Exiting auto-correlation generator ..." << std::endl;
            }
            else {
                // Pass through pixels, assign pixels to buckets 
                if (_verbose) std::cout << "Entering cross-correlation generator ..." << std::endl;
                // The key is a global bucketgrid index and the value is a 
                // list of indices that represent points inside that bucket
                BucketToPoints bucketPointsMapA, bucketPointsMapB;
                // The key is a global bucketgrid index and the value is
                // a list of neighboring buckets
                BucketToBuckets bucketNeighborsMap;
                // First pass through the data, fill buckets with point indices.
                // Also create a lookup tables for points->buckets and buckets->neighbors
                std::vector<double> position(3);
                std::vector<int> binNeighbors;
                for(int i = 0; i < _a.size(); ++i) {
                    position[0] = _a.at(i).x;
                    position[1] = _a.at(i).y;
                    position[2] = _a.at(i).z;
                    int bucketIndex = _bucketgrid.getIndex(position);
                    if(bucketPointsMapA.count(bucketIndex) > 0) {
                        bucketPointsMapA[bucketIndex].push_back(i);
                    }
                    else {
                        bucketPointsMapA[bucketIndex] = std::vector<int>(1,i);
                        _bucketgrid.getBinNeighbors(bucketIndex, binNeighbors);
                        bucketNeighborsMap[bucketIndex] = binNeighbors;
                    }
                }

                for(int i = 0; i < _b.size(); ++i) {
                    position[0] = _b.at(i).x;
                    position[1] = _b.at(i).y;
                    position[2] = _b.at(i).z;
                    int bucketIndex = _bucketgrid.getIndex(position);
                    if(bucketPointsMapB.count(bucketIndex) > 0) {
                        bucketPointsMapB[bucketIndex].push_back(i);
                    }
                    else {
                        bucketPointsMapB[bucketIndex] = std::vector<int>(1,i);
                        _bucketgrid.getBinNeighbors(bucketIndex, binNeighbors);
                        bucketNeighborsMap[bucketIndex] = binNeighbors;
                    }
                }

                int nbuckets = bucketNeighborsMap.size();
                std::cout << "We have " << nbuckets << " buckets" << std::endl;

                // Loop over all buckets
                for(auto &bucket : bucketPointsMapA){
                    // Loop over all points in each bucket
                    for(int i : bucket.second) {
                        // Compare this points to all points in neighboring buckets
                        for(int bucketindex : bucketNeighborsMap[bucket.first]) {
                            // Loop over all points in neighboring bucket
                            for(int j : bucketPointsMapB[bucketindex]) {
                                yield(PairType(_a.at(i),_b.at(j)));
                            }
                        }
                    }
                }
                if (_verbose) std::cout << "Exiting cross-correlation generator ..." << std::endl;
            }
        }
    private:
        PixelIterable _a, _b;
        lk::BinnedGrid _bucketgrid;
        bool _verbose;
        bool _auto;
    }; // BucketPairSearch

    typedef BucketPairSearch<Pixels> BucketSearch;

} // turbooctospice

#endif // TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
