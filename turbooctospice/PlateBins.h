// Created 4-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_PLATE_BINS
#define TURBOOCTOSPICE_PLATE_BINS

#include <CCfits/CCfits>

#include "SkyBins.h"

#include "constants.h"
#include "types.h"

#include "boost/lexical_cast.hpp"

namespace turbooctospice {

    /// A map of items that organized by BOSS Plates
    template<typename T> class PlateBins : public SkyBins<T> {
    public:
        typedef std::map<int, SkyObject> PlateList;

        /// Create a new PlateBins object
        PlateBins(std::string platelist_filename) : SkyBins<T>() {
            std::cout << "Plate bin angular radius: "
                << boost::lexical_cast<std::string>(1.5*deg2rad) << std::endl;
            try {
                // Lookup file name for this target
                std::cout << "Reading platelist file " << platelist_filename << std::endl;
                // Read fits file
                std::unique_ptr<CCfits::FITS> platelist_file(new CCfits::FITS(platelist_filename, CCfits::Read));
                // Read header keywords
                CCfits::PHDU& header = platelist_file->pHDU();
                // Read table columns
                CCfits::ExtHDU& table = platelist_file->extension(1);

                std::vector<int> plates;
                std::vector<float> ra_centers, dec_centers;

                table.column("PLATE").read(plates, 1, table.rows());
                table.column("RACEN").read(ra_centers, 1, table.rows());
                table.column("DECCEN").read(dec_centers, 1, table.rows());

                for(int i = 0; i < plates.size(); ++i) {
                    platelist_[plates[i]] = SkyObject(ra_centers[i]*deg2rad, dec_centers[i]*deg2rad);
                }
            }
            catch (CCfits::FitsException& e) {
                std::cerr << "CCfits Exception Thrown :" << e.message();
            }
        };
        virtual ~PlateBins() {};
        /// Add an item's index to the bin containing the specified angular position
        /// @param ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param item The item to add
        void addItem(const Forest &forest, const T &item) {
            // Find the healpix bin for this ra,dec and save it's index
            int binIndex(forest.plate);

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
            SkyObject pointing(ra, dec);
            for(const auto &plate : platelist_) {
                double cos_sep = pointing.angularSeparation(plate.second);
                double separation = std::acos(cos_sep);
                if(std::fabs(separation - 1.5*deg2rad) < radius){
                    neighbors.push_back(plate.first);
                }
            }
            return neighbors;
        };
        void printBinCenter(const int index) const {
            SkyObject p = platelist_.at(index);
            std::cout << boost::lexical_cast<std::string>(p.dec) << " " << boost::lexical_cast<std::string>(p.ra) << std::endl;
        }
    private:
        PlateList platelist_;
    }; // PlateBins
    typedef PlateBins<int> PlateBinsI;

} // turbooctospice

#endif // TURBOOCTOSPICE_PLATE_BINS
