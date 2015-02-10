// Created 4-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HEALPIX_BINS
#define TURBOOCTOSPICE_HEALPIX_BINS

#include "healpix_map.h"

namespace turbooctospice {

    template<typename T> class HealpixBins {
    public:
        typedef std::map<int, std::vector<T> > Bins;
        typedef Healpix_Map<double> HealpixMap;
        HealpixBins(int order) {
            _map = HealpixMap(order, RING);
        };
        // add an item's index to the bin containing the specified angular position
        void addItem(double theta, double phi, const T &item) {
            // Find the healpix bin for this quasar and save it's index
            int binIndex(_map.ang2pix( {theta, phi} ));
            if(_bins.count(binIndex) > 0) {
                _bins[binIndex].push_back(item);
            }
            else {
                _bins[binIndex] = std::vector<T>(1, item);
            }
        };
        // Return bin indices within specified radius of specified angular position
        std::vector<int> getBinIndicesWithinRadius(double theta, double phi, double radius) const {
            rangeset<int> r;
            _map.query_disc_inclusive({theta, phi}, radius, r);
            std::vector<int> neighbors;
            r.toVector(neighbors);
            return neighbors;
        };
        // Return true if the specified bin index has any contents
        bool checkBinExists(int index) const {
            if(_bins.find(index) == _bins.end()) return false;
            return true;
        };
        // Return the contents of a bin
        const typename Bins::mapped_type& getBinContents(const typename Bins::key_type& k) const {
            return _bins.at(k);
        };
        // Return the Healpix bin index containing the specified angular position
        int ang2pix(double theta, double phi) const { return _map.ang2pix( {theta, phi} ); };
        // Return the number of bins that contain at least one item
        int getNBins() const { return _bins.size(); };
    private:
        int _order;
        Bins _bins;
        HealpixMap _map;
    }; // HealpixBins

    typedef HealpixBins<int> HealpixBinsI;

} // turbooctospice

#endif