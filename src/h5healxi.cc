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

    int order, nthreads, num_sightlines;
    double OmegaLambda, OmegaMatter;
    std::string infile, outfile, axis1, axis2, axis3, covfile, platelist_filename;
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
            "Axis-3 binning, log10(1+z)")
        ("polar", "(r,mu,z) binning")
        ("cart", "(r_perp,r_par,z) binning")
        ("debug", "debug flag")
        ("save-subsamples", "save xi estimates on individual subsamples")
        ("num-sightlines,n", po::value<int>(&num_sightlines)->default_value(0),
            "number of sightlines to use.")
        ("platelist-filename", po::value<std::string>(&platelist_filename)->default_value("/data/sas/dr12/boss/spectro/redux/v5_7_0/platelist.fits"),
            "platelist filename")
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
    bool verbose(vm.count("verbose")), polar(vm.count("polar")), cart(vm.count("cart")),
        debug(vm.count("debug")), save_subsamples(vm.count("save-subsamples"));

    if (infile.size() == 0) {
        std::cerr << "Must specify input file." << std::endl;
        return -1;
    }

    // set up cosmology
    cosmo::AbsHomogeneousUniversePtr cosmology;
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    try {
        cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda, OmegaMatter));
    }
    catch(cosmo::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    // create grid for binning
    tos::AbsTwoPointGridPtr grid;
    tos::XiEstimator::BinningCoordinateType type;

    try {
        lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2), bins3 = lk::createBinning(axis3);
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
            grid.reset(new tos::QuasarGrid(bins1, bins2, bins3));
        }
        if(verbose) {
            std::cout << "Binning grid configuration: " << std::endl;
            for(int axis = 0; axis < 3; ++axis) {
                std::cout << " axis" << axis+1 << " : " << grid->getAxisMin(axis) << " " << grid->getAxisMax(axis) << " "
                    << grid->getAxisBinWidth(axis) << " " << grid->getAxisNBins(axis) << std::endl;
            }
        }
    }
    catch(likely::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    try {
        // load forest sight lines
        tos::HDF5Delta file(infile);
        if (verbose) {
            std::cout << "Loading forest delta fields from " << infile << std::endl;
        }
        std::vector<tos::Forest> sightlines(file.loadForests());
        // trim sample, if requested
        if(num_sightlines > 0 && num_sightlines < sightlines.size()){
            std::vector<tos::Forest>::const_iterator first = sightlines.begin();
            std::vector<tos::Forest>::const_iterator last = sightlines.begin() + num_sightlines;
            std::vector<tos::Forest> selection(first, last);
            sightlines = selection;
        }
        num_sightlines = sightlines.size();

        // create sky bins
        tos::SkyBinsIPtr skybins;
        // which sky binning are we using?
        if(order > 0) {
            skybins.reset(new tos::HealpixBinsI(order));
        }
        else {
            skybins.reset(new tos::PlateBinsI(platelist_filename));
        }

        // add sightlines to sky bins
        unsigned num_pixels(0);
        double min_loglam(sightlines[0].pixels[0].loglam), max_loglam(sightlines[0].pixels[0].loglam);
        for(auto &sightline : sightlines) {
            int npixels = sightline.pixels.size();
            num_pixels += npixels;
            // add sightline to sky bins
            skybins->addItem(sightline, sightline.forest_id);
            // find minimum loglam
            auto loglam_first(sightline.pixels[0].loglam);
            if(loglam_first < min_loglam) {
                min_loglam = loglam_first;
            }
            // find maximum loglam
            auto loglam_last(sightline.pixels[npixels-1].loglam);
            if(loglam_last > max_loglam) {
                max_loglam = loglam_last;
            }
        }

        // sightline/pixel stats
        std::cout << "Read " << num_pixels << " from " << num_sightlines << " lines of sight (LOS)" << std::endl;
        if(verbose) {
            double avg_pixels_per_los(static_cast<double>(num_pixels)/num_sightlines);
            std::cout << "Average number of pixels per LOS: " << boost::lexical_cast<std::string>(avg_pixels_per_los) << std::endl;
            unsigned num_skybins_occupied(skybins->getNBins());
            std::cout << "Number of sky bins occupied: " << num_skybins_occupied << std::endl;
        }
        // the minimum redshift sets the angular scale we will need to consider
        double zmin(std::pow(10, min_loglam-tos::logLyA)-1);
        double max_transverse_scale(cosmology->getTransverseComovingScale(zmin));
        if(verbose) {
            double max_ang(grid->maxAngularScale(max_transverse_scale));
            std::cout << "Transverse comoving scale at z = " << boost::lexical_cast<std::string>(zmin) <<  " : "
                << boost::lexical_cast<std::string>(max_transverse_scale) << " (Mpc/h)" << std::endl;
            std::cout << "Max angular scale at z = " << boost::lexical_cast<std::string>(zmin) <<  " : "
                << boost::lexical_cast<std::string>(max_ang) << " (rad)" << std::endl;
            // for curiosity's sake
            double zmax(std::pow(10, max_loglam-tos::logLyA)-1);
            double min_transverse_scale(cosmology->getTransverseComovingScale(zmax));
            double min_ang(grid->maxAngularScale(min_transverse_scale));
            std::cout << "Transverse comoving scale at z = " << boost::lexical_cast<std::string>(zmax) <<  " : "
                << boost::lexical_cast<std::string>(min_transverse_scale) << " (Mpc/h)" << std::endl;
            std::cout << "Max angular scale at z = " << boost::lexical_cast<std::string>(zmax) <<  " : "
                << boost::lexical_cast<std::string>(min_ang) << " (rad)" << std::endl;
            std::cout << std::endl;
        }

        // run the estimator
        tos::XiEstimator xiest(max_transverse_scale, grid, type, sightlines, skybins);
        xiest.run(nthreads);
        std::cout << std::endl;
        xiest.save_results(outfile);
        if(save_subsamples) {
            xiest.save_subsamples(outfile);
        }
    }
    catch(tos::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}
