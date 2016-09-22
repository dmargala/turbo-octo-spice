// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_TYPES
#define TURBOOCTOSPICE_TYPES

#include "boost/smart_ptr.hpp"
#include "RuntimeError.h"

#include <cmath>
#include <vector>

namespace turbooctospice {

    /// A simple representation of an object in the sky.
    struct SkyObject {
        double ra, dec, sin_dec, cos_dec, sin_ra, cos_ra;
        /// Create a new SkyObject
        /// @param ra_ Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param dec_ Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        SkyObject(double ra_, double dec_) :
        ra(ra_), dec(dec_) {
            sin_dec = std::sin(dec);
            cos_dec = std::cos(dec);
            sin_ra = std::sin(ra);
            cos_ra = std::cos(ra);
        }
        /// Copy constructor
        /// @param other The other SkyObject
        SkyObject(const SkyObject& other) {
            ra = other.ra; dec = other.dec;
            sin_dec = other.sin_dec; cos_dec = other.cos_dec;
            sin_ra = other.sin_ra; cos_ra = other.cos_ra;
        };
        SkyObject() {};
        /// Return cosine of the angular separation between this and the provided SkyObject
        /// @param other The other SkyObject
        double angularSeparation(SkyObject const &other) const {
            return sin_dec*other.sin_dec + cos_dec*other.cos_dec*(sin_ra*other.sin_ra + cos_ra*other.cos_ra);
        }
    };

    /// Represents a LyA forest pixel
    struct ForestPixel {
        float value, loglam, weight, distance;
        /// Create a ForestPixel
        /// @param val The "value" for this pixel
        /// @param loglambda The log10 wavelength of this pixel
        /// @param wgt The weight for this pixel
        /// @param dist The comoving distance to this pixel
        ForestPixel(float val, float loglambda, float wgt, float dist)
        : value(val), loglam(loglambda), weight(wgt), distance(dist) {};
    };

    /// Represents a LyA forest sight line
    struct Forest : SkyObject {
        std::vector<ForestPixel> pixels;
        int forest_id, plate;
        /// Create a Forest sightline
        /// @param _ra Right ascension. The angular distance of a point east of the First Point of Aries, measured along the celestial equator, in radians (0 < ra < 2*pi).
        /// @param _dec Declination. The angular distance of a point north or south of the celestial equator, in radians (-pi/2 < dec < pi/2).
        /// @param _forest_id A unique identifier for this forest. Often used internally and assigned incrementally (0->num_forests).
        /// @param _plate The BOSS plate id containing this sightline.
        Forest(double _ra, double _dec, int _forest_id, int _plate) :
        SkyObject(_ra, _dec), forest_id(_forest_id), plate(_plate) {
        };
    };

    /// Represents a correlation function bin
    struct XiBin {
        double didj, di, dj, wgt, wi, wj, z;
        long num_pairs;
        /// Initialize XiBin
        XiBin() : didj(0), di(0), dj(0), wgt(0), wi(0), wj(0), z(0), num_pairs(0) {};
        /// Accumulate a ForestPixel pair
        /// @param i First ForestPixel
        /// @param j Second ForestPixel
        void accumulate_pair(ForestPixel const &i, ForestPixel const &j) {
            float wdi(i.weight*i.value), wdj(j.weight*j.value);
            didj += wdi*wdj;
            di += wdi;
            dj += wdj;
            wgt += i.weight*j.weight;
            wi += i.weight;
            wj += j.weight;
            z += 0.5*(i.loglam+j.loglam);
            ++num_pairs;
        }
        /// Return the weighted mean value of this XiBin.
        double getMeanProduct() {
            if(wgt <= 0) {
                throw RuntimeError("XiBin: wgt <= 0");
            }
            return didj/wgt;
        }
        /// Finalize the XiBin by normalizing sums.
        void finalize() {
            (wgt > 0) ? didj /= wgt : didj = 0;
            (wi > 0) ? di /= wi : di = 0;
            (wj > 0) ? dj /= wj : dj = 0;
            (num_pairs > 0) ? z /= num_pairs : z = 0;
        }
        /// Add rhs to this XiBin
        /// @param rhs XiBin
        XiBin& operator+=(const XiBin& rhs) {
            didj += rhs.didj;
            di += rhs.di;
            dj += rhs.dj;
            wgt += rhs.wgt;
            wi += rhs.wi;
            wj += rhs.wj;
            z += rhs.z;
            num_pairs += rhs.num_pairs;
            return *this;
        }
        /// Adds two XiBins
        /// @param lhs XiBin
        /// @param rhs XiBin
        friend XiBin operator+(XiBin lhs, const XiBin& rhs){
            return lhs += rhs;
        }
    };

    class AbsTwoPointGrid;
    typedef boost::shared_ptr<AbsTwoPointGrid> AbsTwoPointGridPtr;

    template<typename T> class SkyBins;
    typedef boost::shared_ptr<SkyBins<int> > SkyBinsIPtr;
}

#endif
