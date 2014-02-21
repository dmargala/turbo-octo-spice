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

const double lyA = 1216;

int main(int argc, char **argv) {

    int order;
    double OmegaLambda, OmegaMatter, zmin, zmax, rmax, combine, speclo, foresthi, forestlo;
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
        ("forest-lo", po::value<double>(&forestlo)->default_value(1050),
            "Lyman-alpha forest low cutoff wavelength")
        ("forest-hi", po::value<double>(&foresthi)->default_value(1200),
            "Lyman-alpha forest high cutoff wavelength")
        ("spec-lo", po::value<double>(&speclo)->default_value(3650),
            "Spectrograph wavelength lower limit")
        ("z-max", po::value<double>(&zmax)->default_value(3.5),
            "Maximum redshift to consider")
        ("combine", po::value<double>(&combine)->default_value(4),
            "Number of wavelength bins to combine in fake spectra.")

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
    std::vector<double> sth,cth,sph,cph,s;
    std::vector<std::vector<double> > spectra;

    int goodz = 0;

    double m = std::pow(std::pow(10,0.0001),combine);

    long sumlength = 0;
    long totalpixels = 0;

    for(int i = 0; i < columns[0].size(); ++i) {
        double ra(deg2rad*columns[0][i]);
        double dec(deg2rad*columns[1][i]);
        double z = columns[2][i];
        if(z < zmin || z > zmax) continue;

        double zlo = std::max(speclo/lyA-1,forestlo/lyA*(1+z)-1);
        double zhi = foresthi/lyA*(1+z)-1;

        double zpix = zlo;
        std::vector<double> spectrum;
        while(zpix < zhi) {
            spectrum.push_back(cosmology->getLineOfSightComovingDistance(zpix));
            zpix = (1 + zpix)*m - 1;
            totalpixels++;
        }
        spectra.push_back(spectrum);

        sumlength += spectrum.size();

        double theta = (90.0*deg2rad-dec);

        s.push_back(cosmology->getLineOfSightComovingDistance(z));
        sth.push_back(std::sin(theta));
        cth.push_back(std::cos(theta));
        sph.push_back(std::sin(ra));
        cph.push_back(std::cos(ra));

        pointing p(theta, ra);
        pointings.push_back(p);
        int index = map.ang2pix(p);
        if(buckets.count(index) > 0) {
            buckets[index].push_back(goodz++);
        }
        else {
            buckets[index] = std::vector<int>(1,goodz++);
        }
    }

    std::cout << "Average forest size: " <<  float(sumlength)/spectra.size() <<  " pixels" << std::endl;

    int nbuckets = buckets.size();
    std::cout << "We have " << nbuckets << " buckets w/ data" << std::endl;

    long npair(0), nused(0), npixelpairs(0), npixelpairsused(0);

    rangeset<int> neighbors_rangeset;
    std::vector<int> neighbors;

    double rmaxsq = rmax*rmax;
    double rminsq = 0;
    double cosmax = std::cos(maxAng);

    BOOST_FOREACH(BucketToPixels::value_type &pixels, buckets) {
        // Loop over all points in each bucket
        BOOST_FOREACH(int i, pixels.second) {
            double s_i = s[i];
            double sth_i = sth[i];
            double cth_i = cth[i];
            double sph_i = sph[i];
            double cph_i = cph[i];
            double s_i_sq = s_i*s_i;

            map.query_disc_inclusive(pointings[i], maxAng, neighbors_rangeset);
            neighbors_rangeset.toVector(neighbors);
            // Compare this point to all points in neighboring buckets
            BOOST_FOREACH(int neighbor, neighbors) {
                // Loop over all points in neighboring bucket
                BOOST_FOREACH(int j, buckets[neighbor]) {
                    // Only count pairs once
                    if(j <= i) continue;
                    npair++;
                    double cosij = sth_i*sth[j]*(cph_i*cph[j] + sph_i*sph[j]) + cth_i*cth[j];
                    if(cosij < cosmax) continue;
                    BOOST_FOREACH(double ipix, spectra[i]){
                        BOOST_FOREACH(double jpix, spectra[j]){
                            npixelpairs++;
                            double distsq = ipix*ipix + jpix*jpix - 2*ipix*jpix*cosij;
                            if(distsq >= rmaxsq) continue;
                            if(distsq < rminsq) continue;
                            npixelpairsused++;
                        }
                    }
                    nused++;
                }
            }
        }
    }

    long n(pointings.size());
    long ndistinct = (n*(n-1))/2;
    double consideredFrac = float(npair)/ndistinct;
    double usedFrac = float(nused)/npair;
    std::cout << "Number of distinct los pairs " << ndistinct << std::endl;
    std::cout << "considered " << npair << " of distinct los pairs. (" << consideredFrac << ")" << std::endl;
    std::cout << "used " << nused << " of los pairs considered. (" << usedFrac << ")" << std::endl;

    long ndistinctpixels = (totalpixels*(totalpixels-1))/2;
    double consideredPixelsFrac = float(npixelpairs)/ndistinctpixels;
    double usedPixelsFrac = float(npixelpairsused)/npixelpairs;
    std::cout << "Number of distinct pixel pairs " << ndistinctpixels << std::endl;
    std::cout << "considered " << npixelpairs << " of distinct pixel pairs. (" << consideredPixelsFrac << ")" << std::endl;
    std::cout << "used " << npixelpairsused << " of pixel pairs considered. (" << usedPixelsFrac << ")" << std::endl;

    return 0;
}