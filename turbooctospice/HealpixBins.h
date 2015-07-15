// Created 4-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HEALPIX_BINS
#define TURBOOCTOSPICE_HEALPIX_BINS

#include "healpix_map.h"

#include "constants.h"
#include "types.h"

#include "SkyBins.h"

#include <map>

namespace turbooctospice {

    /// A map of items that organized by HEALPixels
    template<typename T> class HealpixBins : public SkyBins<T> {
    public:
        typedef Healpix_Map<double> HealpixMap;
        /// Create a new HealpixBins object
        /// @param order The HEALPix resolution order
        HealpixBins(int order) : SkyBins<T>() {
            map_ = HealpixMap(order, RING);
            std::cout << "Healpix bin angular \"radius\": " << map_.max_pixrad() << std::endl;
        };
        /// Add an item's index to the bin containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param item The item to add
        void addItem(const Forest &forest, const T &item) {
            // Find the healpix bin for this quasar and save it's index
            int binIndex(ang2pix(forest.ra, forest.dec));

            if(SkyBins<T>::bins_.count(binIndex) > 0) {
                SkyBins<T>::bins_[binIndex].push_back(item);
            }
            else {
                SkyBins<T>::bins_[binIndex] = std::vector<T>(1, item);
                SkyBins<T>::occupied_bins_.push_back(binIndex);
            }
            ++SkyBins<T>::num_entries_;
        };
        /// Return bin indices within radius of an angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param radius The angular radius to query in radians.
        /// @param fact HEALPix search precision parameter
        std::vector<int> getBinIndicesWithinRadius(double ra, double dec, double radius) const {
            std::vector<int> neighbors;
            int fact(4);
            map_.query_disc_inclusive(radec2pnt(ra, dec), radius, neighbors, fact);
            return neighbors;
        };
        /// Return the Healpix bin index containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        int ang2pix(double ra, double dec) const {
            try {
                return map_.ang2pix(radec2pnt(ra, dec));
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
        /// Print bin center of specified index
        void printBinCenter(const int index) const {
            auto p = map_.pix2ang(index);
            std::cout << 0.5*pi - p.theta << " " << p.phi;
        }
    private:
        int order_;
        HealpixMap map_;
    }; // HealpixBins

    typedef HealpixBins<int> HealpixBinsI;

} // turbooctospice

#endif // TURBOOCTOSPICE_HEALPIX_BINS
