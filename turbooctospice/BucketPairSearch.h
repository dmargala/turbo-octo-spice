// Created 04-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
#define TURBOOCTOSPICE_BUCKET_PAIR_SEARCH

#include "likely/likely.h"

#include "boost/foreach.hpp"

#include <iostream>


namespace lk = likely;

namespace turbooctospice {
    template <typename PixelIterable>
    class BucketPairSearch : private PixelIterable {
    public:
        BucketPairSearch(lk::BinnedGrid const &bucketgrid) {
            _bucketgrid = bucketgrid;
        }
        ~BucketPairSearch() {};
        void findPairs(PairGenerator::caller_type& yield, PixelIterable const &a, PixelIterable const &b) const {
            // Pass through pixels, assign pixels to buckets 

            // loop over buckets 
                // Loop over all points in each bucket
                    // Compare this points to all points in neighboring buckets
                        // Loop over all points in neighboring bucket

        }
        void findPairs(PairGenerator::caller_type& yield, PixelIterable const &a) const {
            // The key is a global bucketgrid index and the value is a 
            // list of indices that represent points inside that bucket
            typedef std::map<int, std::vector<int> > BucketIndexToIntegerList;
            BucketIndexToIntegerList bucketPointsMap;
            // The key is a global bucketgrid index and the value is
            // a list of neighboring buckets
            BucketIndexToIntegerList bucketNeighborsMap;
            // First pass through the data, fill buckets with point indices.
            // Also create a lookup tables for points->buckets and buckets->neighbors
            std::vector<double> position(3);
            std::vector<int> binNeighbors;
            for(int i = 0; i < a.size(); ++i) {
                position[0] = a[i].x;
                position[1] = a[i].y;
                position[2] = a[i].z;
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

            int nbuckets = bucketPointsMap.size();
            std::cout << "We have " << nbuckets << " buckets" << std::endl;

            // Loop over all buckets
            BOOST_FOREACH(BucketIndexToIntegerList::value_type &bucket, bucketPointsMap){
                // Loop over all points in each bucket
                BOOST_FOREACH(int i, bucket.second) {
                    // Compare this points to all points in neighboring buckets
                    BOOST_FOREACH(int b, bucketNeighborsMap[bucket.first]) {
                        // Loop over all points in neighboring bucket
                        BOOST_FOREACH(int j, bucketPointsMap[b]) {
                            // Only count pairs once
                            if(j <= i) continue;
                            yield(a[i],a[j]);
                        }
                    }
                }
            }
        }
    private:
        lk::BinnedGrid _bucketgrid;
    }; // BucketPairSearch
} // turbooctospice

#endif // TURBOOCTOSPICE_BUCKET_PAIR_SEARCH
