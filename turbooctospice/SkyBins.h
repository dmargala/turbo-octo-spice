// Created 4-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_SKY_BINS
#define TURBOOCTOSPICE_SKY_BINS

#include <map>

#include "constants.h"
#include "types.h"

#include "boost/lexical_cast.hpp"

namespace turbooctospice {

    /// A map of items that organized by Sky Regions
    template<typename T> class SkyBins {
    public:
        typedef std::map<int, std::vector<T> > Bins;

        /// Create a new SkyBins object
        SkyBins() :
        num_entries_(0) {
            bins_ = Bins();
            occupied_bins_ = std::vector<int>();
        };
        virtual ~SkyBins() {};
        /// Add an item's index to the bin containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param item The item to add
        virtual void addItem(const Forest &forest, const T &item) = 0;
        std::vector<int> getOccupiedBins() const {
            return occupied_bins_;
        }
        /// Return bin indices within radius of an angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param radius The angular radius to query in radians.
        /// @param fact HEALPix search precision parameter
        virtual std::vector<int> getBinIndicesWithinRadius(double ra, double dec, double radius) const = 0;
        /// Return true if the specified bin index has any contents
        /// @param index PlateBins index
        bool checkBinExists(int index) const {
            if(bins_.find(index) == bins_.end()) return false;
            return true;
        };
        /// Return the contents of a bin
        /// @param k Map key
        const typename Bins::mapped_type& getBinContents(const typename Bins::key_type& k) const {
            return bins_.at(k);
        };
        /// Return the number of bins that contain at least one item
        int getNBins() const { return bins_.size(); };
        /// Return the total number of entries.
        unsigned long getNEntries() const { return num_entries_; };
        /// Print bin center of specified index
        virtual void printBinCenter(const int index) const = 0;
    protected:
        Bins bins_;
        unsigned long num_entries_;
        std::vector<int> occupied_bins_;
    private:
    }; // SkyBins

    typedef SkyBins<int> SkyBinsI;

} // turbooctospice

#endif // TURBOOCTOSPICE_SKY_BINS
