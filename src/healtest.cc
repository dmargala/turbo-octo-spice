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
    double OmegaLambda, OmegaMatter, zmin, rmax;
    std::string infile;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("order", po::value<int>(&order)->default_value(1),
            "Healpix map order parameter")
        ("nest", "Use nest ordering scheme for Healpix map")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ("z-min", po::value<double>(&zmin)->default_value(2.1),
            "Minimum z value, sets spherical bin surface distance")
        ("r-max", po::value<double>(&rmax)->default_value(200),
            "Maximum r value to bin.")
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
    bool verbose(vm.count("verbose")), useNest(vm.count("nest"));

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

    double scale = cosmology->getTransverseComovingScale(zmin);
    std::cout << "Transverse comoving scale at z = 2.1: " << scale << std::endl;
    double maxAng = rmax/scale;
    std::cout << "Maximum distance at z = 2.1 (rad): " << maxAng << std::endl;

    Healpix_Ordering_Scheme scheme = (useNest ? NEST : RING);
    Healpix_Map<double> map(order, scheme); 

    std::cout << "Number of pixels: " << map.Npix() << std::endl;
    std::cout << "Max ang dist between any pixel center and its corners: \n\t" 
        << map.max_pixrad() << " rad (" << map.max_pixrad()*scale << " Mpc/h)" << std::endl;

    const double pi = std::atan(1.0)*4;
    const double deg2rad = pi/180.;

    typedef std::map<int, std::vector<int> > BucketToPixels;
    BucketToPixels buckets;
    std::vector<pointing> pointings;
    std::vector<double> X, Y, Z;

    for(int i = 0; i < columns[0].size(); ++i) {
        double ra(deg2rad*columns[0][i]), dec(deg2rad*columns[1][i]);
        double theta = (90.0*deg2rad-dec);
        pointing p(theta, ra);
        pointings.push_back(p);


        double s(cosmology->getLineOfSightComovingDistance(std::abs(columns[2][i])));
        double cosDEC(std::cos(dec));
        X.push_back(s*cosDEC*std::cos(ra));
        Y.push_back(s*cosDEC*std::sin(ra));
        Z.push_back(s*std::sin(dec));

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

    long npair(0), nused(0);

    rangeset<int> neighbors_rangeset;
    std::vector<int> neighbors;

    double maxRange = 1 - std::cos(maxAng);
    double rmaxsq = rmax*rmax;

    BOOST_FOREACH(BucketToPixels::value_type &pixels, buckets) {
        // Loop over all points in each bucket
        BOOST_FOREACH(int i, pixels.second) {
            double xi = X[i];
            double yi = Y[i];
            double zi = Z[i];
            map.query_disc_inclusive(pointings[i], maxAng, neighbors_rangeset);
            neighbors_rangeset.toVector(neighbors);
            // Compare this point to all points in neighboring buckets
            BOOST_FOREACH(int neighbor, neighbors) {
                // Loop over all points in neighboring bucket
                BOOST_FOREACH(int j, buckets[neighbor]) {
                    // Only count pairs once
                    if(j <= i) continue;
                    npair++;
                    double dx = xi - X[j];
                    double dy = yi - Y[j];
                    double dz = zi - Z[j];
                    double dist = dx*dx + dy*dy + dz*dz;
                    if(dist >= rmaxsq) continue;
                    if(dist < 0) continue;
                    nused++;
                }
            }
        }
    }

    long n(columns[0].size());
    long ndistinct = (n*(n-1))/2;
    double consideredFrac = npair*1.0/ndistinct;
    double usedFrac = nused*1.0/npair;
    std::cout << "Number of distinct pairs " << ndistinct << std::endl;
    std::cout << "considered " << npair << " of distinct pairs. (" << consideredFrac << ")" << std::endl;
    std::cout << "used " << nused << " of pairs considered. (" << usedFrac << ")" << std::endl;

    return 0;
}