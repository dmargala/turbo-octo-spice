// Created 20-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Healpix test

#include "healpix_map.h"

#include "boost/program_options.hpp"
#include "boost/foreach.hpp"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include <iostream>
#include <fstream>
#include <map>
#include <vector>

namespace po = boost::program_options;
namespace lk = likely;

int main(int argc, char **argv) {

    int order;
    double OmegaLambda, OmegaMatter;
    std::string infile;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("order", po::value<int>(&order)->default_value(1),
            "Healpix map order parameter")
        ("ring", "Use ring ordering scheme for Healpix map")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
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
    bool verbose(vm.count("verbose")), useRing(vm.count("ring"));

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

    double unitBAO = 200/scale;
    std::cout << "Approx BAO scale at z = 2.1 in rad: " << unitBAO << std::endl;

    Healpix_Ordering_Scheme scheme = (useRing ? RING : NEST);
    Healpix_Map<double> map(order, scheme); 

    std::cout << "Number of pixels: " << map.Npix() << std::endl;

    std::cout << "Max ang dist (in radian) between any pixel center and its corners: \n\t" << map.max_pixrad() << std::endl;

    const double PI = std::atan(1.0)*4;
    const double Deg2Rad = PI/180.;

    typedef std::map<int, std::vector<int> > BucketToPixels;
    BucketToPixels buckets;
    std::vector<pointing> pointings;
    for(int i = 0; i < columns[0].size(); ++i) {
        double phi = columns[0][i]*Deg2Rad;
        double theta = (90.-columns[1][i])*Deg2Rad;
        pointing p(theta, phi);
        pointings.push_back(p);

        int index = map.ang2pix(p);
        if(buckets.count(index) > 0) {
            buckets[index].push_back(i);
        }
        else {
            buckets[index] = std::vector<int>(1,i);
        }
    }

    int nbuckets = buckets.size();
    std::cout << "We have " << nbuckets << " buckets w/ data" << std::endl;

    long npair = 0;

    rangeset<int> neighbors_rangeset;
    std::vector<int> neighbors;

    BOOST_FOREACH(BucketToPixels::value_type &pixels, buckets) {
        // Loop over all points in each bucket
        BOOST_FOREACH(int i, pixels.second) {
            rangeset<int> neighbors_rangeset;
            map.query_disc(pointings[i], unitBAO, neighbors_rangeset);
            neighbors_rangeset.toVector(neighbors);
            // Compare this points to all points in neighboring buckets
            BOOST_FOREACH(int neighbor, neighbors) {
                // Loop over all points in neighboring bucket
                BOOST_FOREACH(int j, buckets[neighbor]) {
                    if(j == i) std::cout << "We found ourself!" << std::endl;
                    // Only count pairs once
                    if(j <= i) continue;
                    npair++;
                }
            }
        }
    }

    std::cout << "Number of pairs counted: " << npair++ << std::endl;
    // fix_arr<int, 8> neighbors;

    // map.neighbors(0, neighbors);

    // for(int i = 0; i < neighbors.size(); ++i) {
    //  std::cout << neighbors[i] << std::endl;
    // }

    return 0;
}