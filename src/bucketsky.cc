// Created 30-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Bucket correlation function estimator.

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"
#include "boost/foreach.hpp"

#include "boost/geometry.hpp"
#include "boost/geometry/geometries/point.hpp"
#include "boost/geometry/geometries/box.hpp"

#include "boost/geometry/index/rtree.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

#include "tos.h"

namespace po = boost::program_options;
namespace lk = likely;

namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

namespace tos = turbooctospice;


//typedef bg::model::point<float, 3, bg::cs::cartesian> point;

// class Triangle {

// public:
//     Triangle(int v1, int v2, int v3) : 
//     _v1(v1), _v2(v2), _v3(v3) {

//     }
//     int _v1, _v2, _v3;

// private:

// };

// const double t = (1.0 + std::sqrt(5.0)) / 2.0;

// class TriangleMesh {

// public:

//     std::vector<tos::Triangle> triangles;
//     std::vector<point> points;
//     std::map<long, int> midPointCache;
//     int index;
    
//     TriangleMesh(int recursionLevel){

//         index = 0;

//         // create 12 vertices of a icosahedron
//         addVertex(point(-1, t, 0));
//         addVertex(point( 1, t, 0));
//         addVertex(point(-1,-t, 0));
//         addVertex(point( 1, t, 0));

//         addVertex(point( 0,-1, t));
//         addVertex(point( 0, 1, t));
//         addVertex(point( 0,-1,-t));
//         addVertex(point( 0, 1,-t));

//         addVertex(point( t, 0,-1));
//         addVertex(point( t, 0, 1));
//         addVertex(point(-t, 0,-1));
//         addVertex(point(-t, 0, 1));

//         // create 20 triangles of the icosahedron
//         // 5 faces around point 0
//         triangles.push_back(tos::Triangle(0, 11, 5));
//         triangles.push_back(tos::Triangle(0, 5, 1));
//         triangles.push_back(tos::Triangle(0, 1, 7));
//         triangles.push_back(tos::Triangle(0, 7, 10));
//         triangles.push_back(tos::Triangle(0, 10, 11));

//         // 5 adjacent faces 
//         triangles.push_back(tos::Triangle(1, 5, 9));
//         triangles.push_back(tos::Triangle(5, 11, 4));
//         triangles.push_back(tos::Triangle(11, 10, 2));
//         triangles.push_back(tos::Triangle(10, 7, 6));
//         triangles.push_back(tos::Triangle(7, 1, 8));

//         // 5 faces around point 3
//         triangles.push_back(tos::Triangle(3, 9, 4));
//         triangles.push_back(tos::Triangle(3, 4, 2));
//         triangles.push_back(tos::Triangle(3, 2, 6));
//         triangles.push_back(tos::Triangle(3, 6, 8));
//         triangles.push_back(tos::Triangle(3, 8, 9));

//         // 5 adjacent faces 
//         triangles.push_back(tos::Triangle(4, 9, 5));
//         triangles.push_back(tos::Triangle(2, 4, 11));
//         triangles.push_back(tos::Triangle(6, 2, 10));
//         triangles.push_back(tos::Triangle(8, 6, 7));
//         triangles.push_back(tos::Triangle(9, 8, 1));

//         // refine triangles
//         for (int i = 0; i < recursionLevel; ++i) {
//             std::vector<tos::Triangle> newTriangles;
//             BOOST_FOREACH(tos::Triangle tri, triangles) {
//                 // replace triangle by 4 triangles
//                 int a = addMidPoint(tri._v1, tri._v2);
//                 int b = addMidPoint(tri._v2, tri._v3);
//                 int c = addMidPoint(tri._v3, tri._v1);

//                 newTriangles.push_back(tos::Triangle(tri._v1, a, c));
//                 newTriangles.push_back(tos::Triangle(tri._v2, b, a));
//                 newTriangles.push_back(tos::Triangle(tri._v3, c, b));
//                 newTriangles.push_back(tos::Triangle(a, b, c));
//             }
//             triangles = newTriangles;
//         }

//     };

//     int addVertex(point p) {
//         double norm = bg::distance(p,point(0,0,0));
//         points.push_back(point(bg::get<0>(p)/norm,bg::get<1>(p)/norm,bg::get<2>(p)/norm));
//         return index++;
//     }

//     int addMidPoint(int p1, int p2) {
//         bool firstIsSmaller = p1 < p2;
//         long smallerIndex = firstIsSmaller ? p1 : p2;
//         long greaterIndex = firstIsSmaller ? p2 : p1;
//         long key = (smallerIndex << 32) + greaterIndex;

//         if(midPointCache.count(key) > 0){
//             return midPointCache[key];
//         }

//         point point1(points[p1]);
//         point point2(points[p2]);
//         point midpoint(
//             bg::get<0>(point1)+bg::get<0>(point2),
//             bg::get<1>(point1)+bg::get<1>(point2),
//             bg::get<2>(point1)+bg::get<2>(point2));
//         int i = addVertex(midpoint);
//         midPointCache[key] = i;
//         return i;
//     };

// private:

// };


int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string infile,outfile;
    double OmegaLambda,OmegaMatter;
    int recursionLevel;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ("output,o", po::value<std::string>(&outfile)->default_value(""),
            "Filename to write to")
        ("level", po::value<int>(&recursionLevel)->default_value(0),
            "recursion level to use for triangulation")
        ;

    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > columns(3);
    try {
        std::ifstream in(infile.c_str());
        lk::readVectors(in,columns);
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << infile << ": " << e.what() << std::endl;
        return -3;
    }
    if(verbose) {
        std::cout << "Read " << columns[0].size() << " rows from " << infile
            << std::endl;
    }

    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmo::AbsHomogeneousUniversePtr cosmology(
        new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));

    double scale = cosmology->getTransverseComovingScale(2.1);
    std::cout << "Transverse comoving scale at z = 2.1: " << scale << std::endl;

    tos::TriangleMesh tmesh(recursionLevel);
    int ntriangles = tmesh.triangles.size();
    int npoints = tmesh.points.size();
    std::cout << "Triangular mesh has " << ntriangles << " triangles and " << npoints << " vertices." << std::endl;

    double sidelength = bg::distance(tmesh.points[npoints-1], tmesh.points[npoints-2]);
    std::cout << "Approx triangle side length: " << sidelength << std::endl;

    typedef bg::model::point<double, 2, bg::cs::spherical_equatorial<bg::degree> > spherical_point;
    typedef bg::model::box<spherical_point> box;
    typedef std::pair<spherical_point, int> value;

    // dec = acos(z/r)
    // ra = atan2(y,x)

    // create the rtree using default constructor
    typedef bgi::rtree< value, bgi::rstar<16> > Rtree;
    Rtree rtree;
    for(int i = 0; i < columns[0].size(); ++i) {
        spherical_point p(columns[0][i], columns[1][i]);
        box b = bg::return_envelope<box>(p);
        rtree.insert(std::make_pair(p, i));
        //std::cout << "Distance : " << distance(p0, p)*scale << std::endl;
    }


    double unitBAO = 110./scale;
    std::cout << "Approx BAO scale at z = 2.1 in rad: " << unitBAO << std::endl;

    double limit(3*unitBAO);

    long count = 0;
    for(int i = 0; i < columns[0].size()-1; ++i) {
        double ra = columns[0][i];
        double dec = columns[1][i];
        spherical_point p(ra,dec);
        box query_box(spherical_point(ra - unitBAO, dec - limit), spherical_point(ra+limit, dec+limit));

        std::vector<value> result_s;
        rtree.query(bgi::within(query_box), std::back_inserter(result_s));

        // BOOST_FOREACH(value const&v, result_s) {
        //     int j = v.second;
        //     if (j <= i) continue;
        //     spherical_point q(v.first);
        //     double dist = scale*bg::distance(p,q);
        //     if (dist > 200) continue;
        //     count++;
        // }
    }

    std::cout << count << std::endl;

    // // Save the estimator results
    // try {
    //     std::ofstream out(outfile.c_str());
    //     for(int index = 0; index < xi.size(); ++index) {
    //         out << index << ' ' << xi[index] << std::endl;
    //     }
    //     out.close();
    // }
    // catch(std::exception const &e) {
    //     std::cerr << "Error while saving results: " << e.what() << std::endl;
    // }

    return 0;
}
