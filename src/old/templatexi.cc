// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Example strategy (polcicy) design pattern using templated classes

#include "boost/program_options.hpp"

#include "turbooctospice.h"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include <iostream>
#include <fstream>
#include <vector>

namespace po = boost::program_options;
namespace lk = likely;
namespace tos = turbooctospice;

const double lyA = 1216;

int main(int argc, char **argv) {

    // Configure command-line option processing
    double OmegaLambda, OmegaMatter, zmin, zmax, combine, speclo, foresthi, forestlo;
    std::string infile,outfile,axis,buckets;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value(""),
            "Filename to write correlation function to")
        ("axis", po::value<std::string>(&axis)->default_value("[0:200]*50"),
            "Xi axis binning")
        ("ignore", "No binnig (use to profile pair search methods)")
        ("buckets", po::value<std::string>(&buckets)->default_value(""),
            "Specify bucket search binning scheme")
        ("norm", "Normalize xi by dividing by weights")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.728,"0.728"),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("z-min", po::value<double>(&zmin)->default_value(2.1,"2.1"),
            "Minimum z value, sets spherical bin surface distance")
        ("z-max", po::value<double>(&zmax)->default_value(3.5),
            "Maximum redshift to consider")
        ("forest-lo", po::value<double>(&forestlo)->default_value(1040),
            "Lyman-alpha forest low cutoff wavelength")
        ("forest-hi", po::value<double>(&foresthi)->default_value(1200),
            "Lyman-alpha forest high cutoff wavelength")
        ("spec-lo", po::value<double>(&speclo)->default_value(3650),
            "Spectrograph wavelength lower limit")
        ("combine", po::value<double>(&combine)->default_value(4),
            "Number of wavelength bins to combine in fake spectra.")
        ("heal", "Use HEALPIX search")
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
    bool verbose(vm.count("verbose")),ignore(vm.count("ignore")),norm(vm.count("norm")),heal(vm.count("heal"));

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

    // Instantiate the correlation function grid
    lk::AbsBinningCPtr bins = lk::createBinning(axis);
    int nbins(bins->getNBins());
    double min(bins->getBinLowEdge(0)), max(bins->getBinHighEdge(nbins-1));

    // Run the estimator
    std::vector<double> xi;

    if (heal) {
        if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
        cosmo::AbsHomogeneousUniversePtr cosmology(
            new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));

        double scale = cosmology->getTransverseComovingScale(zmin);
        double maxAng = max/scale;
        if(verbose) {
            std::cout << "Transverse comoving scale at z = 2.1: " << scale << std::endl;
            std::cout << "Maximum distance at z = 2.1 (rad): " << maxAng << std::endl;
        }

        // Healpix_Ordering_Scheme scheme = RING
        int order = 5;
        Healpix_Map<double> map(order, RING); 

        if(verbose) {
            std::cout << "Number of pixels: " << map.Npix() << std::endl;
            std::cout << "Max ang dist between any pixel center and its corners: \n\t" 
                << map.max_pixrad() << " rad (" << map.max_pixrad()*scale << " Mpc/h)" << std::endl;
        }

        const double pi = std::atan(1.0)*4;
        const double deg2rad = pi/180.;

        long totalpixels = 0;

        tos::Quasars quasars;
        double m = std::pow(std::pow(10,0.0001),combine);

        for(int i = 0; i < columns[0].size(); ++i){
            double ra(deg2rad*columns[0][i]);
            double dec(deg2rad*columns[1][i]);
            double theta = (90.0*deg2rad-dec);
            double z = columns[2][i];
            if(z < zmin || z > zmax) continue;

            double zlo = std::max(speclo/lyA-1,forestlo/lyA*(1+z)-1);
            double zhi = foresthi/lyA*(1+z)-1;

            double zpix = zlo;
            tos::Quasarf quasar;
            quasar.p = pointing(theta, ra);
            float sth = std::sin(theta);
            float cth = std::cos(theta);
            float sph = std::sin(ra);
            float cph = std::cos(ra);
            float s;
            std::vector<tos::LOSPixelf> pixels;
            while(zpix < zhi) {
                s = cosmology->getLineOfSightComovingDistance(zpix);
                pixels.push_back(tos::LOSPixelf(s,sth,cth,sph,cph,1,1));
                zpix = (1 + zpix)*m - 1;
            }
            totalpixels += pixels.size();
            quasar.pixels = pixels;
            quasars.push_back(quasar);
        }

        long nquasars = quasars.size();
        long ndistinct = nquasars*(nquasars-1)/2;
        long ndistinctpixels = totalpixels*(totalpixels-1)/2;

        std::cout << "Average forest size: " <<  ((double)totalpixels)/nquasars <<  " pixels" << std::endl;
        std::cout << "Number of distinct los pairs " << ndistinct << std::endl;
        std::cout << "Number of distinct pixel pairs " << ndistinctpixels << std::endl;

        if(ignore) {
            tos::HealIgnoreXi xiestimator(
                new tos::HealSearch(quasars, map, maxAng, verbose), 
                new tos::IgnoreAng, 
                verbose);
            xi = xiestimator.run(norm);
        }
        else {
            tos::HealBinXi xiestimator(
                new tos::HealSearch(quasars, map, maxAng, verbose), 
                new tos::BinAng(min, max, nbins), 
                verbose);
            xi = xiestimator.run(norm);
        }

    }
    else {
        // Parse the input file
        tos::Pixels pixels;
        for(int i = 0; i < columns[0].size(); ++i) {
            tos::Pixel pixel;
            pixel.x = columns[0][i];
            pixel.y = columns[1][i];
            pixel.z = columns[2][i];
            pixel.d = 1;//columns[3][i];
            pixel.w = 1;//columns[4][i];
            pixel.i = i;
            pixels.push_back(pixel);
        }
        if(buckets.length()) {
            lk::AbsBinningCPtr bucketbins1 = lk::createBinning(buckets),
                bucketbins2 = lk::createBinning(buckets),
                bucketbins3 = lk::createBinning(buckets);
            lk::BinnedGrid bucketgrid(bucketbins1, bucketbins2, bucketbins3);

            if(ignore) {
                tos::BucketIgnoreXi xiestimator(
                    new tos::BucketSearch(pixels, bucketgrid, verbose), 
                    new tos::Ignore, 
                    verbose);
                xi = xiestimator.run(norm); 
            }
            else {
                tos::BucketBinXi xiestimator(
                    new tos::BucketSearch(pixels, bucketgrid,  verbose), 
                    new tos::BinXYZ(min, max, nbins), 
                    verbose);
                xi = xiestimator.run(norm); 
            }
        }
        else {
            if(ignore) {
                tos::BruteIgnoreXi xiestimator(
                    new tos::BruteSearch(pixels, verbose), 
                    new tos::Ignore, 
                    verbose);
                xi = xiestimator.run(norm); 
            }
            else {
                tos::BruteBinXi xiestimator(
                    new tos::BruteSearch(pixels, verbose), 
                    new tos::BinXYZ(min, max, nbins), 
                    verbose);
                xi = xiestimator.run(norm); 
            }
        }
    }

    // Save the estimator results
    if(outfile.length() > 0) {
        if(verbose) {
            std::cout << "Saving xi to " << outfile << std::endl;
        }
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
    }
}