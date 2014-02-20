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

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string infile,outfile;
    double OmegaLambda,OmegaMatter;
    int recursionOrder;
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
        ("order", po::value<int>(&recursionOrder)->default_value(0),
            "recursion order to use for triangulation")
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

    tos::TriangleMesh tmesh(recursionOrder);
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
