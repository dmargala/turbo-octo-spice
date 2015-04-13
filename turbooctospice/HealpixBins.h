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
            std::cout << "Max pix rad: " << _map.max_pixrad() << std::endl;
        };
        // add an item's index to the bin containing the specified angular position
        void addItem(double theta, double phi, const T &item) {
            // Find the healpix bin for this quasar and save it's index
            int binIndex(ang2pix(theta, phi));

            if(_bins.count(binIndex) > 0) {
                _bins[binIndex].push_back(item);
            }
            else {
                _bins[binIndex] = std::vector<T>(1, item);
            }
            ++_nentries;
        };
        // Return bin indices within radius of an angular position
        std::vector<int> getBinIndicesWithinRadius(double theta, double phi, double radius, int fact = 4) const {
            std::vector<int> neighbors;
            _map.query_disc_inclusive({theta, phi}, radius, neighbors, fact);
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
        int ang2pix(double theta, double phi) const { 
            try {
                return _map.ang2pix( {theta, phi} ); 
            }
            catch(PlanckError e) {
                std::cerr << "Invalid (theta,phi): " << theta << ", " << phi << std::endl;
                return 0;
            }
        };
        // Return the number of bins that contain at least one item
        int getNBins() const { return _bins.size(); };
        unsigned long getNEntries() const { return _nentries; };
        void printBinCenter(const int index) const {
            auto p = _map.pix2ang(index);
            std::cout << p.theta << " " << p.phi;
        }
    private:
        int _order;
        Bins _bins;
        HealpixMap _map;
        unsigned long _nentries;
    }; // HealpixBins

    typedef HealpixBins<int> HealpixBinsI;

} // turbooctospice

#endif