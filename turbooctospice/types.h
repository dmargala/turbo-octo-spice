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

    	// Pixel(float _x, float _y, float _z, float _d, float _w) : x(_x), y(_y), z(_z), d(_d), w(_w) {}
    };

    typedef std::vector<Pixel> Pixels;

    typedef std::pair<Pixel,Pixel> PixelPair;
    typedef boost::coroutines::coroutine<PixelPair()> PairGenerator;
}

#endif