// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_TRIANGLE_MESH
#define TURBOOCTOSPICE_TRIANGLE_MESH

#include "types.h"
#include "Triangle.h"

namespace turbooctospice {

	class TriangleMesh {

	public:
	    std::vector<Triangle> triangles;
	    std::vector<point> points;
	    std::map<long, int> midPointCache;
	    int index;
	    
	    TriangleMesh(int recursionLevel);
	    int addVertex(point p);
	    int addMidPoint(int p1, int p2);

	private:

	};

}

#endif
