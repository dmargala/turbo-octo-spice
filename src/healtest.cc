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

typedef std::map<int, std::vector<int> > BucketToPixels;

typedef double dtype;

struct QuasarSpectrum {
    dtype s, sth, cth, sph, cph;
    pointing p;
    std::vector<dtype> spectrum;
    std::vector<dtype> lambdas;
    std::vector<dtype> weights;
};

void healxi(Healpix_Map<double> const &map, BucketToPixels &buckets, std::vector<QuasarSpectrum> &quasars,
lk::BinnedGrid const &grid, bool rmu, double maxAng, double x1min, double x1max, double x2min, double x2max, std::vector<double> &xi) {
    // create internal accumulation vectors
    double x1minSq(x1min*x1min), x1maxSq(x1max*x1max);
    double cosmax = std::cos(maxAng);

    int nbins = grid.getNBinsTotal();
    std::vector<double> dsum(nbins,0.), wsum(nbins,0.);

    rangeset<int> neighbors_rangeset;
    std::vector<int> neighbors;

    long nuniquepairs(0), nlospairsused(0), totalpixels(0), npixelpairs(0), npixelpairsused(0);

    long neighborBucketCounter(0), neighborQuasarCounter(0);

    for(int i = 0; i < quasars.size(); ++i) {
        QuasarSpectrum qi = quasars[i];

        totalpixels += qi.spectrum.size();

        map.query_disc_inclusive(qi.p, maxAng, neighbors_rangeset);
        neighbors_rangeset.toVector(neighbors);
        // Compare this point to all points in neighboring buckets
        BOOST_FOREACH(int neighbor, neighbors) {
            neighborBucketCounter++;
            // Loop over all points in neighboring bucket
            BOOST_FOREACH(int j, buckets[neighbor]) {
                neighborQuasarCounter++;
                // Only count pairs once
                if(j <= i) continue;
                nuniquepairs++;
                QuasarSpectrum qj = quasars[j];
                double cosij = qi.sth*qj.sth*(qi.cph*qj.cph + qi.sph*qj.sph) + qi.cth*qj.cth;
                if(cosij < cosmax) continue;
                nlospairsused++;
                for(int ipix = 0; ipix < qi.spectrum.size(); ++ipix) {
                    dtype si = qi.lambdas[ipix];
                    dtype di = qi.spectrum[ipix];
                    dtype wi = qi.weights[ipix];
                    for(int jpix = 0; jpix < quasars[j].spectrum.size(); ++jpix) {
                        npixelpairs++;
                        dtype sj = qj.lambdas[jpix];
                        double distSq = si*si + sj*sj - 2*si*sj*cosij;
                        if(distSq < x1minSq || distSq > x1maxSq) continue;
                        npixelpairsused++;
                    }
                }
            }
        }
    }

    // Write quasar pair statistics to console
    long nquasars(quasars.size());
    long ndistinct = (nquasars*(nquasars-1))/2;
    double consideredFrac = float(nuniquepairs)/ndistinct;
    double usedFrac = float(nlospairsused)/nuniquepairs;
    std::cout << "neighboring buckets: " << neighborBucketCounter << std::endl;
    std::cout << "neighboring quasars: " << neighborQuasarCounter << std::endl;
    std::cout << "Number of distinct los pairs " << ndistinct << std::endl;
    std::cout << "considered " << nuniquepairs << " of distinct los pairs. (" << consideredFrac << ")" << std::endl;
    std::cout << "used " << nlospairsused << " of los pairs considered. (" << usedFrac << ")" << std::endl;
    // Write pixel pair statistics to console
    long ndistinctpixels = (totalpixels*(totalpixels-1))/2;
    double consideredPixelsFrac = float(npixelpairs)/ndistinctpixels;
    double usedPixelsFrac = float(npixelpairsused)/npixelpairs;
    std::cout << "Number of distinct pixel pairs " << ndistinctpixels << std::endl;
    std::cout << "considered " << npixelpairs << " of distinct pixel pairs. (" << consideredPixelsFrac << ")" << std::endl;
    std::cout << "used " << npixelpairsused << " of pixel pairs considered. (" << usedPixelsFrac << ")" << std::endl;

    for(int index = 0; index < nbins; ++index) {
        if(wsum[index] > 0) dsum[index] /= wsum[index];
    }
    dsum.swap(xi);
}

int main(int argc, char **argv) {

    int order;
    double OmegaLambda, OmegaMatter, zmin, zmax, rmax, combine, speclo, foresthi, forestlo;
    std::string infile, outfile, axis1, axis2;
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
        ("forest-lo", po::value<double>(&forestlo)->default_value(1040),
            "Lyman-alpha forest low cutoff wavelength")
        ("forest-hi", po::value<double>(&foresthi)->default_value(1200),
            "Lyman-alpha forest high cutoff wavelength")
        ("spec-lo", po::value<double>(&speclo)->default_value(3650),
            "Spectrograph wavelength lower limit")
        ("z-max", po::value<double>(&zmax)->default_value(3.5),
            "Maximum redshift to consider")
        ("combine", po::value<double>(&combine)->default_value(4),
            "Number of wavelength bins to combine in fake spectra.")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
            "Axis-1 binning")
        ("axis2", po::value<std::string>(&axis2)->default_value("[0:200]*50"),
            "Axis-2 binning")
        ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")

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
    bool verbose(vm.count("verbose")), useNest(vm.count("nest")), rmu(vm.count("rmu"));

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

    BucketToPixels buckets;
    std::vector<QuasarSpectrum> quasars;

    int goodz = 0;

    double m = std::pow(std::pow(10,0.0001),combine);

    long totalpixels = 0;

    for(int i = 0; i < columns[0].size(); ++i) {
        double ra(deg2rad*columns[0][i]);
        double dec(deg2rad*columns[1][i]);
        double theta = (90.0*deg2rad-dec);
        double z = columns[2][i];
        if(z < zmin || z > zmax) continue;

        double zlo = std::max(speclo/lyA-1,forestlo/lyA*(1+z)-1);
        double zhi = foresthi/lyA*(1+z)-1;

        double zpix = zlo;
        std::vector<dtype> spectrum;
        std::vector<dtype> weights;
        std::vector<dtype> lambdas;
        while(zpix < zhi) {
            lambdas.push_back(cosmology->getLineOfSightComovingDistance(zpix));
            weights.push_back(1);
            spectrum.push_back(1);
            zpix = (1 + zpix)*m - 1;
        }
        totalpixels += spectrum.size();

        QuasarSpectrum qso;
        qso.spectrum = spectrum;
        qso.weights = weights;
        qso.lambdas = lambdas;
        qso.s = cosmology->getLineOfSightComovingDistance(z);
        qso.sth = std::sin(theta);
        qso.cth = std::cos(theta);
        qso.sph = std::sin(ra);
        qso.cph = std::cos(ra);
        //pointing p(theta, ra);
        qso.p = pointing(theta, ra);

        quasars.push_back(qso);
        int index = map.ang2pix(qso.p);
        if(buckets.count(index) > 0) {
            buckets[index].push_back(goodz++);
        }
        else {
            buckets[index] = std::vector<int>(1,goodz++);
        }
    }

    std::cout << "Average forest size: " <<  double(totalpixels)/quasars.size() <<  " pixels" << std::endl;

    int nbuckets = buckets.size();
    std::cout << "We have " << nbuckets << " buckets w/ data" << std::endl;

    // Generate the correlation function grid and run the estimator
    std::vector<double> xi;
    try {
        lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2);
        double x1min(bins1->getBinLowEdge(0)), x1max(bins1->getBinHighEdge(bins1->getNBins()-1));
        double x2min(bins2->getBinLowEdge(0)), x2max(bins2->getBinHighEdge(bins2->getNBins()-1));
        lk::BinnedGrid grid(bins1,bins2);
        healxi(map, buckets,quasars,grid,rmu,maxAng,x1min,x1max,x2min,x2max,xi);
    }
    catch(std::exception const &e) {
        std::cerr << "Error while running the estimator: " << e.what() << std::endl;
    }

    // Save the estimator results
    try {
        std::ofstream out(outfile.c_str());
        for(int index = 0; index < xi.size(); ++index) {
            out << index << ' ' << xi[index] << std::endl;
        }
        out.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while saving results: " << e.what() << std::endl;
    }


    return 0;
}