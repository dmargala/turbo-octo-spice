// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_TYPES
#define TURBOOCTOSPICE_TYPES

// #include "boost/geometry.hpp"
// #include "boost/geometry/geometries/point.hpp"

// #include "boost/coroutine/coroutine.hpp"

#include "boost/smart_ptr.hpp"
#include "HealpixBins.h"
#include "RuntimeError.h"

#include <cmath>
#include <vector>

namespace turbooctospice {
	// Represents a 3D point in Cartesisna coordinates.
	// typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;

	// struct QuasarPixel {
 //    	float frac, lam, wgt, dist;
	// };

 //    class Pixel {
 //    public:
 //    	float x, y, z, d, w;
 //    	int i;

 //        /// WWHYYYYY???
 //        bool operator==(const Pixel &rhs) const {
 //            return true;
 //        }
 //        bool operator!=(const Pixel &rhs) const {
 //            return !this->operator==(rhs);
 //        }
 //    };

 //    typedef std::vector<Pixel> Pixels;

     struct SkyObject {
         double ra, dec, sin_dec, cos_dec, sin_ra, cos_ra;
         SkyObject(double ra_, double dec_) :
         ra(ra_), dec(dec_) {
             sin_dec = std::sin(dec);
             cos_dec = std::cos(dec);
             sin_ra = std::sin(ra);
             cos_ra = std::cos(ra);
         }
         SkyObject(const SkyObject& other) {
             ra = other.ra; dec = other.dec;
             sin_dec = other.sin_dec; cos_dec = other.cos_dec;
             sin_ra = other.sin_ra; cos_ra = other.cos_ra;
         };
         SkyObject() {};

         double angularSeparation(SkyObject const &other) const {
             return sin_dec*other.sin_dec + cos_dec*other.cos_dec*(sin_ra*other.sin_ra + cos_ra*other.cos_ra);
         }
     };

    /// Represents a LyA forest pixel
    struct ForestPixel {
        float value, loglam, weight, distance;
        ForestPixel(float val, float loglambda, float wgt, float dist)
        : value(val), loglam(loglambda), weight(wgt), distance(dist) {};
    };

    /// Represents a LyA forest sight line
    struct Forest : SkyObject {
        std::vector<ForestPixel> pixels;
        int forest_id, plate;
        Forest(double _ra, double _dec, int _forest_id, int _plate) :
        SkyObject(_ra, _dec), forest_id(_forest_id), plate(_plate) {
        };
    };

    /// Represents a correlation function bin
    struct XiBin {
        double didj, di, dj, wgt, wi, wj;
        long num_pairs;
        XiBin() : didj(0), di(0), dj(0), wgt(0), wi(0), wj(0), num_pairs(0) {};
        void accumulate_pair(ForestPixel const &i, ForestPixel const &j) {
            float wdi(i.weight*i.value), wdj(j.weight*j.value);
            didj += wdi*wdj;
            di += wdi;
            dj += wdj;
            wgt += i.weight*j.weight;
            wi += i.weight;
            wj += j.weight;
            ++num_pairs;
        }
        double getMeanProduct() {
            if(wgt <= 0) {
                throw RuntimeError("XiBin: wgt <= 0");
            }
            return didj/wgt;
        }
        void finalize() {
            (wgt > 0) ? didj /= wgt : didj = 0;
            (wi > 0) ? di /= wi : di = 0;
            (wj > 0) ? dj /= wj : dj = 0;
        }
        XiBin& operator+=(const XiBin& rhs) {
            didj += rhs.didj;
            di += rhs.di;
            dj += rhs.dj;
            wgt += rhs.wgt;
            wi += rhs.wi;
            wj += rhs.wj;
            num_pairs += rhs.num_pairs;
            return *this;
        }
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
