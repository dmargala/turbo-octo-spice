#include "boost/foreach.hpp"

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point.hpp"

#include "TriangleMesh.h"

#include <cmath>

namespace bg = boost::geometry;

namespace local = turbooctospice;

    
local::TriangleMesh::TriangleMesh(int recursionLevel) {
    index = 0;

    const double goldenRatio((1.0+std::sqrt(5.0))/2.0);

    // create 12 vertices of a icosahedron
    local::TriangleMesh::addVertex(point(-1, goldenRatio, 0));
    addVertex(point( 1, goldenRatio, 0));
    addVertex(point(-1,-goldenRatio, 0));
    addVertex(point( 1, goldenRatio, 0));

    addVertex(point( 0,-1, goldenRatio));
    addVertex(point( 0, 1, goldenRatio));
    addVertex(point( 0,-1,-goldenRatio));
    addVertex(point( 0, 1,-goldenRatio));

    addVertex(point( goldenRatio, 0,-1));
    addVertex(point( goldenRatio, 0, 1));
    addVertex(point(-goldenRatio, 0,-1));
    addVertex(point(-goldenRatio, 0, 1));

    // create 20 triangles of the icosahedron
    // 5 faces around point 0
    triangles.push_back(local::Triangle(0, 11, 5));
    triangles.push_back(local::Triangle(0, 5, 1));
    triangles.push_back(local::Triangle(0, 1, 7));
    triangles.push_back(local::Triangle(0, 7, 10));
    triangles.push_back(local::Triangle(0, 10, 11));

    // 5 adjacent faces 
    triangles.push_back(local::Triangle(1, 5, 9));
    triangles.push_back(local::Triangle(5, 11, 4));
    triangles.push_back(local::Triangle(11, 10, 2));
    triangles.push_back(local::Triangle(10, 7, 6));
    triangles.push_back(local::Triangle(7, 1, 8));

    // 5 faces around point 3
    triangles.push_back(local::Triangle(3, 9, 4));
    triangles.push_back(local::Triangle(3, 4, 2));
    triangles.push_back(local::Triangle(3, 2, 6));
    triangles.push_back(local::Triangle(3, 6, 8));
    triangles.push_back(local::Triangle(3, 8, 9));

    // 5 adjacent faces 
    triangles.push_back(local::Triangle(4, 9, 5));
    triangles.push_back(local::Triangle(2, 4, 11));
    triangles.push_back(local::Triangle(6, 2, 10));
    triangles.push_back(local::Triangle(8, 6, 7));
    triangles.push_back(local::Triangle(9, 8, 1));

    // refine triangles
    for (int i = 0; i < recursionLevel; ++i) {
        std::vector<local::Triangle> newTriangles;
        BOOST_FOREACH(local::Triangle tri, triangles) {
            // replace triangle by 4 triangles
            int a = local::TriangleMesh::addMidPoint(tri._v1, tri._v2);
            int b = addMidPoint(tri._v2, tri._v3);
            int c = addMidPoint(tri._v3, tri._v1);

            newTriangles.push_back(local::Triangle(tri._v1, a, c));
            newTriangles.push_back(local::Triangle(tri._v2, b, a));
            newTriangles.push_back(local::Triangle(tri._v3, c, b));
            newTriangles.push_back(local::Triangle(a, b, c));
        }
        triangles = newTriangles;
    }

}

int local::TriangleMesh::addVertex(point p) {
    double norm = bg::distance(p,point(0,0,0));
    points.push_back(point(bg::get<0>(p)/norm,bg::get<1>(p)/norm,bg::get<2>(p)/norm));
    return index++;
}

int local::TriangleMesh::addMidPoint(int p1, int p2) {
    bool firstIsSmaller = p1 < p2;
    long smallerIndex = firstIsSmaller ? p1 : p2;
    long greaterIndex = firstIsSmaller ? p2 : p1;
    long key = (smallerIndex << 32) + greaterIndex;

    if(midPointCache.count(key) > 0){
        return midPointCache[key];
    }

    point point1(points[p1]);
    point point2(points[p2]);
    point midpoint(
        bg::get<0>(point1)+bg::get<0>(point2),
        bg::get<1>(point1)+bg::get<1>(point2),
        bg::get<2>(point1)+bg::get<2>(point2));
    int i = addVertex(midpoint);
    midPointCache[key] = i;
    return i;
}