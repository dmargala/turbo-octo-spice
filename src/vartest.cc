// Created 8-Oct-2012 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Stacks many Gaussian random fields on the field maximum (or minimum).

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

namespace po = boost::program_options;
namespace lk = likely;
 
int main(int argc, char **argv) {
    // Configure command-line option processing
    double spacing, xlos, ylos, zlos, binsize, rmin;
    long npairs;
    int nx,ny,nz,seed,nfields,nbins;
    std::string loadPowerFile, prefix;
    po::options_description cli("Stacks many Gaussian random fields on the field maximum (or minimum).");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("spacing", po::value<double>(&spacing)->default_value(2),
            "Grid spacing in Mpc/h.")
        ("nx", po::value<int>(&nx)->default_value(128),
            "Grid size along x-axis.")
        ("ny", po::value<int>(&ny)->default_value(0),
            "Grid size along y-axis (or zero for ny=nx).")
        ("nz", po::value<int>(&nz)->default_value(0),
            "Grid size along z-axis (or zero for nz=ny).")
        ("load-power", po::value<std::string>(&loadPowerFile)->default_value(""),
            "Reads k,P(k) values (in h/Mpc units) to interpolate from the specified filename.")
        ("seed", po::value<int>(&seed)->default_value(510),
            "Random seed to use for GRF.")
        ("nfields", po::value<int>(&nfields)->default_value(1000),
            "Number of fields to stack.")
        ("test","Use the test fft generator.")
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

    // Fill in any missing grid dimensions.
    if(0 == ny) ny = nx;
    if(0 == nz) nz = ny;
    
    double pi(4*std::atan(1));

    if(verbose) {
        std::cout << "Will stack " << nfields << " GRFs with dimensions (x,y,z) = " 
            << boost::format("(%d,%d,%d)") % nx % ny % nz 
            << boost::format(" using %.2f Mpc/h grid spacing.") % spacing << std::endl;
    }

    // Load a tabulated power spectrum for interpolation.
    cosmo::PowerSpectrumPtr power;
    if(0 < loadPowerFile.length()) {
        std::vector<std::vector<double> > columns(2);
        std::ifstream in(loadPowerFile.c_str());
        lk::readVectors(in,columns);
        in.close();
        if(verbose) {
            std::cout << "Read " << columns[0].size() << " rows from " << loadPowerFile
                << std::endl;
        }
        double twopi2(2*pi*pi);
        // rescale to k^3/(2pi^2) P(k)
        for(int row = 0; row < columns[0].size(); ++row) {
            double k(columns[0][row]);
            columns[1][row] *= k*k*k/twopi2;
        }
        // Create an interpolator of this data.
        lk::InterpolatorPtr iptr(new lk::Interpolator(columns[0],columns[1],"cspline"));
        // Use the resulting interpolation function for future power calculations.
        power = lk::createFunctionPtr(iptr);
    }
    else {
        std::cerr << "Missing required load-power filename." << std::endl;
        return -2;
    }

    // Create the generator.
    cosmo::AbsGaussianRandomFieldGeneratorPtr generator;
    if(vm.count("test")) {
        generator.reset(new cosmo::TestFftGaussianRandomFieldGenerator(power, spacing, nx, ny, nz));
    }
    else {
        generator.reset(new cosmo::FftGaussianRandomFieldGenerator(power, spacing, nx, ny, nz));
    }
    if(verbose) {
        std::cout << "Memory size = "
            << boost::format("%.1f Mb") % (generator->getMemorySize()/1048576.) << std::endl;
    }

    // Accumulator
    lk::WeightedAccumulator varStats;

    // Initialize the random number source.
    lk::RandomPtr random = lk::Random::instance();
    random->setSeed(seed);

    for(int ifield = 0; ifield < nfields; ++ifield){
        // Generate Gaussian random field
        generator->generate();
        lk::WeightedAccumulator deltaStats;
        // Fill 1-d, 2-d histograms
        for(int ix = 0; ix < nx; ++ix){
            for(int iy = 0; iy < ny; ++iy){
                for(int iz = 0; iz < nz; ++iz){
                    deltaStats.accumulate(generator->getField(ix,iy,iz));
                }
            }
        }
        std::cout << deltaStats.mean() << " " << deltaStats.variance() << " " << deltaStats.count() << std::endl;
        varStats.accumulate(deltaStats.variance());
    }

    // Print extreme value mean, variance, and count
    std::cout << boost::format("Delta variance mean, variance, counts: %.10g %.10g %d")
        % varStats.mean() % varStats.variance() % varStats.count() << std::endl;

    // Done
    return 0;
}