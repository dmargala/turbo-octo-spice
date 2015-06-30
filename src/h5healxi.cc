// Created 28-Jan-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "boost/program_options.hpp"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "turbooctospice.h"
#include <stdexcept>
// #include <stdio.h>

namespace po = boost::program_options;
namespace lk = likely;
namespace tos = turbooctospice;

void healxi(tos::HealpixBinsI const &healbins, std::vector<tos::Forest> const &sight_lines,
tos::AbsTwoPointGridPtr const &grid, double max_ang, std::vector<tos::XiBin> &xi, std::vector<std::vector<double >> &xi_cov) {

    unsigned long numLOS(sight_lines.size());
    std::cout << "Number of sight lines: " << numLOS << std::endl;

    // stat counters
    unsigned long numHealpixBinsSearched(0), numLOSPairs(0), numLOSPairsUsed(0),
        numPixels(0), numPixelPairs(0), numPixelPairsUsed(0);

    // create internal accumulation vectors
    int nbins = grid->getNBinsTotal();
    std::cout << "Number of bins: " << nbins << std::endl;
    std::vector<double> dsum(nbins,0), wsum(nbins,0);

    // to avoid calling trig functions inside loop
    double cos_max_ang(std::cos(max_ang));

    // allocate temporary vectors before loop
    std::vector<int> neighbors;


    // temporary variables for separation calculation
    float los_pixel_dist_sq, los_pixel_projection_times_two, distSq;
    float weight, product;

    float dist, mu, zpair;

    // bin index variables
    int binIndex, mubin, zbin;

    // axis binning
    float min_dist(grid->getAxisMin(0)), max_dist(grid->getAxisMax(0)), bin_width(grid->getAxisBinWidth(0));
    float min_distsq(min_dist*min_dist), max_distsq(max_dist*max_dist);

    int nbins1(grid->getAxisNBins(1));

    int nbins2(grid->getAxisNBins(2));
    float min_zpair(grid->getAxisMin(2)), max_zpair(grid->getAxisMax(2)), dz(grid->getAxisBinWidth(2));

    boost::progress_display show_progress(sight_lines.size());

    std::map<int, std::vector<tos::XiBin> > healxis;

    for(int i = 0; i < sight_lines.size(); ++i) {
        ++show_progress;
        auto line_of_sight = sight_lines[i];
        numPixels += line_of_sight.pixels.size();
        auto neighbors = healbins.getBinIndicesWithinRadius(line_of_sight.ra, line_of_sight.dec, max_ang);
        int healpix_index = healbins.ang2pix(line_of_sight.ra, line_of_sight.dec);
        if(!(healxis.count(healpix_index) > 0) ) {
            healxis[healpix_index] = std::vector<tos::XiBin>(nbins, {});
        }
        // search neighboring healpix bins
        for(int neighbor : neighbors) {
            ++numHealpixBinsSearched;
            if(!healbins.checkBinExists(neighbor)) continue;
            // look for potential los pairs inside this healpix bins
            for(int j : healbins.getBinContents(neighbor)) {
                // only count pairs once
                if(j <= i) continue;
                ++numLOSPairs;
                auto other_los = sight_lines[j];
                // check angular separation
                double cos_separation(line_of_sight.angularSeparation(other_los));
                if(cos_separation <= cos_max_ang) continue;
                double separation = std::acos(cos_separation);
                ++numLOSPairsUsed;
                // accumulate statistics for pixel pairs
                for(auto& los_pixel : line_of_sight.pixels) {
                    los_pixel_dist_sq = los_pixel.distance*los_pixel.distance;
                    los_pixel_projection_times_two = 2*los_pixel.distance*cos_separation;
                    for(auto& other_pixel : other_los.pixels) {
                        ++numPixelPairs;
                        // check pairs are within our binning grid

                        // check parallel separation
                        distSq = los_pixel_dist_sq + (other_pixel.distance - los_pixel_projection_times_two)*other_pixel.distance;
                        if(distSq >= max_distsq || distSq < min_distsq) continue;
                        dist = std::sqrt(distSq);
                        binIndex = int((dist - min_dist)/bin_width);

                        // todo: check transverse separation
                        mu = (dist == 0 ? 0 : std::fabs(los_pixel.distance-other_pixel.distance)/dist);
                        mubin = int(mu * nbins1);

                        if(mubin < 0) mubin = 0;
                        if(mubin >= nbins1) mubin = nbins1-1;

                        binIndex = mubin + binIndex*nbins1;

                        // todo: check average pair distance
                        zpair = 0.5*(los_pixel.loglam + other_pixel.loglam) - tos::logLyA;
                        zbin = int((zpair-min_zpair)/dz);
                        binIndex = zbin + binIndex*nbins2;

                        if(binIndex < 0 || binIndex >= nbins) {
                            printf("\n");

                            printf("%d %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g\n", binIndex,
                                los_pixel.distance, line_of_sight.ra, line_of_sight.dec,
                                other_pixel.distance, other_los.ra, other_los.dec,
                                cos_separation, separation, dist, mu, zpair);
                            throw std::runtime_error("invalid bin index");
                            continue;
                        }
                        // accumulate pixel pair
                        weight = los_pixel.weight*other_pixel.weight;
                        product = weight*los_pixel.value*other_pixel.value;
                        healxis[healpix_index][binIndex].didj += product;
                        healxis[healpix_index][binIndex].wgt += weight;
                        ++numPixelPairsUsed;
                    }
                }
            }
        }
    }

    // copy data to output vectors
    std::vector<tos::XiBin> xisum(nbins, {});
    for(int index = 0; index < nbins; ++index) {
        for(const auto& pair : healxis) {
            xisum[index].didj += healxis[pair.first][index].didj;
            xisum[index].wgt += healxis[pair.first][index].wgt;
        }
        if(xisum[index].wgt > 0) {
            xisum[index].didj /= xisum[index].wgt;
        }
    }

    // Estimate covariance matrix
    std::vector<std::vector<double> > cov(nbins, std::vector<double>(nbins, 0));
    for(int a = 0; a < nbins; ++a){
        for(int b = 0; b < nbins; ++b) {
            for(const auto& pair : healxis) {
                cov[a][b] += healxis[pair.first][a].didj*healxis[pair.first][b].didj
                    - healxis[pair.first][a].wgt*healxis[pair.first][b].wgt*xisum[a].didj*xisum[b].didj;
            }
            cov[a][b] /= xisum[a].wgt*xisum[b].wgt;
        }
    }

    xisum.swap(xi);
    cov.swap(xi_cov);

    // print stats
    std::cout << "Number of Healpix bins searched: " << numHealpixBinsSearched << std::endl;

    // line of sight pair statistics
    unsigned long numLOSPairsTotal = (numLOS*(numLOS-1))/2;
    double fracLOSPairsConsidered = (double) numLOSPairs / numLOSPairsTotal;
    double fracLOSPairsUsed = (double) numLOSPairsUsed / numLOSPairs;

    std::cout << "Number of distinct los pairs " << numLOSPairsTotal << std::endl;
    std::cout << "considered " << numLOSPairs << " of distinct los pairs. (" << fracLOSPairsConsidered << ")" << std::endl;
    std::cout << "used " << numLOSPairsUsed << " of los pairs considered. (" << fracLOSPairsUsed << ")" << std::endl;

    // pixel pair statistics
    unsigned long numPixelPairsTotal = (numPixels*(numPixels-1))/2;
    double fracPixelPairsConsidered = (double) numPixelPairs / numPixelPairsTotal;
    double fracPixelPairsUsed = (double) numPixelPairsUsed / numPixelPairs;

    std::cout << "Number of distinct pixel pairs " << numPixelPairsTotal << std::endl;
    std::cout << "considered " << numPixelPairs << " of distinct pixel pairs. (" << fracPixelPairsConsidered << ")" << std::endl;
    std::cout << "used " << numPixelPairsUsed << " of pixel pairs considered. (" << fracPixelPairsUsed << ")" << std::endl;
}

int main(int argc, char **argv) {

    int order;
    double OmegaLambda, OmegaMatter;
    std::string infile, outfile, axis1, axis2, axis3, covfile;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("order", po::value<int>(&order)->default_value(5),
            "Healpix map order parameter")
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

    // initialize Healpix bins
    tos::HealpixBinsI healbins(order);

    // load forest sight lines
    std::vector<tos::Forest> sight_lines;

    tos::HDF5Delta file(infile);
    try {
        sight_lines = file.loadForests(!skip_ngc, !skip_sgc);
    }
    catch(tos::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    // add sight lines to healpix bins
    unsigned long totalpixels(0);
    double min_loglam(sight_lines[0].pixels[0].loglam), max_loglam(sight_lines[0].pixels[0].loglam);
    for(int i = 0; i < sight_lines.size(); ++i) {
        auto los = sight_lines[i].pixels;

        totalpixels += los.size();
        healbins.addItem(sight_lines[i].ra, sight_lines[i].dec, i);

        // find minimum loglam
        if(los[0].loglam < min_loglam) {
            min_loglam = los[0].loglam;
        }
        if(los[los.size()-1].loglam > max_loglam) {
            max_loglam = los[los.size()-1].loglam;
        }
    }
    double zmin(std::pow(10, min_loglam-tos::logLyA)-1);
    double zmax(std::pow(10, max_loglam-tos::logLyA)-1);

    int numHealBinsOccupied(healbins.getNBins());
    std::cout << "Read " << totalpixels << " from " << sight_lines.size() << " lines of sight (LOS)" << std::endl;
    std::cout << "Average number of pixels per LOS: " <<  static_cast<double>(totalpixels)/sight_lines.size() << std::endl;
    std::cout << "Number of Healpix bins occupied: " << numHealBinsOccupied
        << " (" << static_cast<double>(numHealBinsOccupied)/(12*std::pow(4, order)) << ")" << std::endl;

    // create grid for binning
    lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2), bins3 = lk::createBinning(axis3);
    tos::AbsTwoPointGridPtr grid;
    if(polar) {
        grid.reset(new tos::PolarGrid(bins1, bins2, bins3));
    }
    else if (cart) {
        std::cerr << "Cartesigna grid implemented yet!" << std::endl;
        return -1;
        // grid.reset(new tos::CartesianGrid(bins1, bins2, bins3));
    }
    else {
        std::cerr << "Observing grid not implemented yet!" << std::endl;
        return -1;
        // grid.reset(new tos::QuasarGrid(bins1, bins2, bins3));
    }

    // the minimum redshift sets the angular scale we will need to consider
    double scale(cosmology->getTransverseComovingScale(zmin));
    double max_ang(grid->maxAngularScale(scale));
    std::cout << "Transverse comoving scale at z = " << zmin <<  " (Mpc/h): " << scale << std::endl;
    std::cout << "Max angular scale at z = " << zmin <<  " (rad): " << max_ang  << std::endl;

    std::cout << "Transverse comoving scale at z = " << zmax <<  " (Mpc/h): " << cosmology->getTransverseComovingScale(zmax) << std::endl;
    std::cout << "Max angular scale at z = " << zmax <<  " (rad): " << grid->maxAngularScale(cosmology->getTransverseComovingScale(zmax))  << std::endl;

    // Generate the correlation function grid and run the estimator
    std::vector<tos::XiBin> xi;
    std::vector<std::vector<double> > cov;
    try {
        healxi(healbins, sight_lines, grid, max_ang, xi, cov);
    }
    catch(std::exception const &e) {
        std::cerr << "Error while running the estimator: " << e.what() << std::endl;
    }

    // Save the estimator results
    try {
        std::vector<double> binCenters(3);

        std::string estimator_filename(outfile + ".dat");
        if(verbose) {
            std::cout << "Saving correlation function to: " << estimator_filename << std::endl;
        }
        std::ofstream estimator_file(estimator_filename.c_str());
        for(int index = 0; index < xi.size(); ++index) {
            grid->getBinCenters(index, binCenters);
            estimator_file << index << ' ' << binCenters[0] << ' ' << binCenters[1] << ' ' << binCenters[2] << ' '
                << xi[index].didj << ' ' << xi[index].di << ' ' << xi[index].dj << ' ' << xi[index].wgt << std::endl;
        }
        estimator_file.close();

        std::string covariance_filename(outfile + ".cov");
        if(verbose) {
            std::cout << "Saving covariance matrix to: " << covariance_filename << std::endl;
        }
        std::ofstream covariance_file(covariance_filename.c_str());
        for(int a = 0; a < cov.size(); ++a) {
            for(int b = 0; b < cov[0].size(); ++b) {
                covariance_file << b + cov[0].size()*a << ' ' << a << ' ' << b << ' ' << cov[a][b] << std::endl;
            }
        }
        covariance_file.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while saving results: " << e.what() << std::endl;
    }

    return 0;
}
