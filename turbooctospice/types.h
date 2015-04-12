// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_TYPES
#define TURBOOCTOSPICE_TYPES

// #include "boost/geometry.hpp"
// #include "boost/geometry/geometries/point.hpp"

// #include "boost/coroutine/coroutine.hpp"

#include "boost/smart_ptr.hpp"

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

    /// Represents a LyA forest pixel 
    struct ForestPixel {
        float value, loglam, weight, distance;
        ForestPixel(float val, float loglambda, float wgt, float dist) 
        : value(val), loglam(loglambda), weight(wgt), distance(dist) {};
    };

    /// Represents a LyA forest sight line
    struct Forest {
        double theta, phi, sdec, cdec, sph, cph;
        std::vector<ForestPixel> pixels;
    };

    class AbsTwoPointGrid;
    typedef boost::shared_ptr<AbsTwoPointGrid> AbsTwoPointGridPtr;


}

#endif