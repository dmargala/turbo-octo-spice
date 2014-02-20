#ifndef TOS_TYPES
#define TOS_TYPES

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point.hpp"

namespace turbooctospice {
	// Represents a 3D point in Cartesisna coordinates.
	typedef boost::geometry::model::point<float, 3, boost::geometry::cs::cartesian> point;
}

#endif