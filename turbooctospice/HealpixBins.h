// Created 4-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HEALPIX_BINS
#define TURBOOCTOSPICE_HEALPIX_BINS

#include "healpix_map.h"

#include "constants.h"
#include <map>

namespace turbooctospice {

    /// A map of items that organized by HEALPixels
    template<typename T> class HealpixBins {
    public:
        typedef std::map<int, std::vector<T> > Bins;
        typedef Healpix_Map<double> HealpixMap;
        /// Create a new HealpixBins object
        /// @param order The HEALPix resolution order
        HealpixBins(int order) {
            _map = HealpixMap(order, RING);
            std::cout << "Max pix rad: " << _map.max_pixrad() << std::endl;
        };
        /// Add an item's index to the bin containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param item The item to add
        void addItem(double ra, double dec, const T &item) {
            // Find the healpix bin for this quasar and save it's index
            int binIndex(ang2pix(ra, dec));

            if(_bins.count(binIndex) > 0) {
                _bins[binIndex].push_back(item);
            }
            else {
                _bins[binIndex] = std::vector<T>(1, item);
                _occupiedBins.push_back(binIndex);
            }
            ++_nentries;
        };
        std::vector<int> getOccupiedBins() const {
            return _occupiedBins;
        }
        /// Return bin indices within radius of an angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param radius The angular radius to query in radians.
        /// @param fact HEALPix search precision parameter
        std::vector<int> getBinIndicesWithinRadius(double ra, double dec, double radius, int fact = 4) const {
            std::vector<int> neighbors;
            _map.query_disc_inclusive(radec2pnt(ra, dec), radius, neighbors, fact);
            return neighbors;
        };
        /// Return true if the specified bin index has any contents
        /// @param index HealpixBins index
        bool checkBinExists(int index) const {
            if(_bins.find(index) == _bins.end()) return false;
            return true;
        };
        /// Return the contents of a bin
        /// @param k Map key
        const typename Bins::mapped_type& getBinContents(const typename Bins::key_type& k) const {
            return _bins.at(k);
        };
        /// Return the Healpix bin index containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        int ang2pix(double ra, double dec) const {
            try {
                return _map.ang2pix(radec2pnt(ra, dec));
            }
            catch(PlanckError e) {
                std::cerr << "HealBins.ang2pix: " << e.what() << std::endl;
                return 0;
            }
        };
        /// Return a HEALPix pointing object corresponding to the specified ra and dec
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        pointing radec2pnt(double ra, double dec) const {
            return {0.5*pi - dec, ra};
        }
        /// Return the number of bins that contain at least one item
        int getNBins() const { return _bins.size(); };
        /// Return the total number of entries.
        unsigned long getNEntries() const { return _nentries; };
        /// Print bin center of specified index
        void printBinCenter(const int index) const {
            auto p = _map.pix2ang(index);
            std::cout << 0.5*pi - p.theta << " " << p.phi;
        }
    private:
        int _order;
        Bins _bins;
        HealpixMap _map;
        unsigned long _nentries;
        std::vector<int> _occupiedBins;
    }; // HealpixBins

    typedef HealpixBins<int> HealpixBinsI;

} // turbooctospice

#endif
