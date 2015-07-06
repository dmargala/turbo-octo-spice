// Created 28-Jan-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include <iostream>
#include <string>
#include <stdexcept>

#include "boost/program_options.hpp"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "turbooctospice.h"

namespace po = boost::program_options;
namespace lk = likely;
namespace tos = turbooctospice;

int main(int argc, char **argv) {

    int order, nthreads;
    double OmegaLambda, OmegaMatter;
    std::string infile, outfile, axis1, axis2, axis3, covfile;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("order", po::value<int>(&order)->default_value(5),
            "Healpix map order parameter")
        ("nthreads", po::value<int>(&nthreads)->default_value(3),
            "Number of threads to use")
        ("omega-lambda", po::value<double>(&OmegaLambda)->default_value(0.73, "0.73"),
            "Present-day value of OmegaLambda.")
        ("omega-matter", po::value<double>(&OmegaMatter)->default_value(0),
            "Present-day value of OmegaMatter or zero for 1-OmegaLambda.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ("output,o", po::value<std::string>(&outfile)->default_value("healxi"),
            "Output filename base")
        //("axis1", po::value<std::string>(&axis1)->default_value("[0:0.05]*25"),
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
            "Axis-1 binning, r_par (Mpc/h), r (Mpc/h), or dloglam")
        //("axis2", po::value<std::string>(&axis2)->default_value("[5:175]*1"),
        ("axis2", po::value<std::string>(&axis2)->default_value("[0:1]*1"),
            "Axis-2 binning, r_perp (Mpc/h), mu (r_par/r), or dtheta (arcmin)")
        ("axis3", po::value<std::string>(&axis3)->default_value("[0.46:.65]*1"),
            "Axis-3 binning, log10(z+1)")
        ("polar", "(r,mu,z) binning")
        ("cart", "(r_perp,r_par,z) binning")
        ("skip-ngc", "only use sgc sight lines")
        ("skip-sgc", "only use ngc sight lines")
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
    bool verbose(vm.count("verbose")), polar(vm.count("polar")), cart(vm.count("cart")), fits(vm.count("fits")),
         skip_ngc(vm.count("skip-ngc")), skip_sgc(vm.count("skip-sgc"));

    if (infile.size() == 0) {
        std::cerr << "Must specify input file." << std::endl;
        return -1;
    }

    if(skip_ngc && skip_sgc) {
        std::cerr << "Can't specify options '--skip-sgc' and '--skip-ngc' together. Use one or neither." << std::endl;
        return -1;
    }

    // set up cosmology
    cosmo::AbsHomogeneousUniversePtr cosmology;
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda, OmegaMatter));

    // create grid for binning
    lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2), bins3 = lk::createBinning(axis3);
    tos::AbsTwoPointGridPtr grid;
    tos::XiEstimator::BinningCoordinateType type;
    if(polar) {
        type = tos::XiEstimator::PolarCoordinates;
        grid.reset(new tos::PolarGrid(bins1, bins2, bins3));
    }
    else if(cart) {
        type = tos::XiEstimator::CartesianCoordinates;
        grid.reset(new tos::CartesianGrid(bins1, bins2, bins3));
    }
    else {
        type = tos::XiEstimator::ObservingCoordinates;
        // std::cerr << "Observing coordinate grid not implemented yet!" << std::endl;
        // return -1;
        grid.reset(new tos::QuasarGrid(bins1, bins2, bins3));
    }

    try {
        tos::XiEstimator xiest(order, infile, cosmology, grid, type, skip_ngc, skip_sgc);
        xiest.run(nthreads);
        xiest.save_results(outfile);
    }
    catch(tos::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }
    catch(std::exception const &e) {
        std::cerr << "Error while running the estimator: " << e.what() << std::endl;
    }

    return 0;
}
