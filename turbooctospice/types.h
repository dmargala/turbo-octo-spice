// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_TYPES
#define TURBOOCTOSPICE_TYPES

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point.hpp"

#include "boost/coroutine/coroutine.hpp"

namespace turbooctospice {
	// Represents a 3D point in Cartesisna coordinates.
	typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;

	struct QuasarPixel {
    	float frac, lam, wgt, dist;
	};

    struct Pixel {
    	float x, y, z, d, w;
    	int i;
    };

    /// This is the Pixels typedef.
    ///
    typedef std::vector<Pixel> Pixels;

}

#endif