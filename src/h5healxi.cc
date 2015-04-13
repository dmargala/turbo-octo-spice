// Created 28-Jan-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include "boost/program_options.hpp"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "turbooctospice.h"

namespace po = boost::program_options;
namespace lk = likely;
namespace tos = turbooctospice;

const double piovertwo = 0.5*tos::pi;

void healxi(tos::HealpixBinsI const &healbins, std::vector<tos::Forest> const &sight_lines,
tos::AbsTwoPointGridPtr const &grid, double max_ang, std::vector<tos::XiBin> &xi) {

    unsigned long numLOS(sight_lines.size());
    std::cout << "Number of sight lines: " << numLOS << std::endl;

    // stat counters
    unsigned long numHealpixBinsSearched(0), numLOSPairs(0), numLOSPairsUsed(0), 
        numPixels(0), numPixelPairs(0), numPixelPairsUsed(0);

    // create internal accumulation vectors
    int nbins = grid->getNBinsTotal();
    std::cout << "Number of bins: " << nbins << std::endl;
    std::vector<tos::XiBin> xisum(nbins, {});
    std::vector<double> dsum(nbins,0), wsum(nbins,0);

    // to avoid calling trig functions inside loop
    double cos_max_ang(std::cos(max_ang));

    // allocate temporary vectors before loop
    std::vector<int> neighbors;

    // temporary variables for separation calculation
    float los_pixel_dist_sq, los_pixel_projection_times_two, distSq, dist, mu;
    float weight, product;

    // bin index variables
    int binIndex;

    // axis 1
    float min_dist(0), max_dist(200), bin_width(4), nbins1(50);
    float min_distsq(min_dist*min_dist), max_distsq(max_dist*max_dist);
    // axis 2
    // double min_mu(0), max_mu(1), dmu(1.0/3.0), nbins2(3);
    // axis 3
    // todo

    boost::progress_display show_progress(sight_lines.size());

    for(int i = 0; i < sight_lines.size(); ++i) {
        ++show_progress;
        auto line_of_sight = sight_lines[i];
        numPixels += line_of_sight.pixels.size();
        auto neighbors = healbins.getBinIndicesWithinRadius(piovertwo-line_of_sight.dec, line_of_sight.ra, max_ang);
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

                        // check transverse separation
                        // mu = std::fabs(pj.distance-pi.distance)/dist;
                        // if(mu >= max_mu || mu <= min_mu) continue;
                        // binIndex = int((mu - min_mu)/dmu) + binIndex*nbins2;

                        // check average pair distance
                        // todo

                        // accumulate pixel pair
                        weight = los_pixel.weight*other_pixel.weight;
                        product = weight*los_pixel.value*other_pixel.value;
                        dsum[binIndex] += product;
                        wsum[binIndex] += weight;
                        ++numPixelPairsUsed;
                    }
                }
            }
        }
    }

    // copy data to output vectors
    for(int index = 0; index < nbins; ++index) {
        if(wsum[index] > 0) {
            xisum[index].didj = dsum[index]/wsum[index];
            xisum[index].wgt = wsum[index];
        }
    }
    xisum.swap(xi);

    // print stats
    std::cout << "Number of Healpix bins searched: " << numHealpixBinsSearched << std::endl;

    // line of sight pair statistics 
    unsigned long numLOSPairsTotal = (numLOS*(numLOS-1))/2;
    double fracLOSPairsConsidered = static_cast<double>(numLOSPairs)/numLOSPairsTotal;
    double fracLOSPairsUsed = static_cast<double>(numLOSPairsUsed)/numLOSPairs;

    std::cout << "Number of distinct los pairs " << numLOSPairsTotal << std::endl;
    std::cout << "considered " << numLOSPairs << " of distinct los pairs. (" << fracLOSPairsConsidered << ")" << std::endl;
    std::cout << "used " << numLOSPairsUsed << " of los pairs considered. (" << fracLOSPairsUsed << ")" << std::endl;

    // pixel pair statistics
    unsigned long numPixelPairsTotal = (numPixels*(numPixels-1))/2;
    double fracPixelPairsConsidered = static_cast<double>(numPixelPairs)/numPixelPairsTotal;
    double fracPixelPairsUsed = static_cast<double>(numPixelPairsUsed)/numPixelPairs;

    std::cout << "Number of distinct pixel pairs " << numPixelPairsTotal << std::endl;
    std::cout << "considered " << numPixelPairs << " of distinct pixel pairs. (" << fracPixelPairsConsidered << ")" << std::endl;
    std::cout << "used " << numPixelPairsUsed << " of pixel pairs considered. (" << fracPixelPairsUsed << ")" << std::endl;
}

int main(int argc, char **argv) {

    int order;
    double OmegaLambda, OmegaMatter;
    std::string infile, outfile, axis1, axis2, axis3;
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
        ("output,o", po::value<std::string>(&outfile)->default_value("healxi_out.txt"),
            "Output filename")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:0.05]*25"),
            "Axis-1 binning, r_par (Mpc/h), r (Mpc/h), or dloglam")
        ("axis2", po::value<std::string>(&axis2)->default_value("[5:175]*1"),
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
    sight_lines = file.loadForests(!skip_ngc, !skip_sgc);

    // add sight lines to healpix bins
    unsigned long totalpixels(0);
    double min_loglam(sight_lines[0].pixels[0].loglam), max_loglam(sight_lines[0].pixels[0].loglam);
    for(int i = 0; i < sight_lines.size(); ++i) {
        auto los = sight_lines[i].pixels;

        totalpixels += los.size();
        healbins.addItem(piovertwo - sight_lines[i].dec, sight_lines[i].ra, i);

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
    lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2), 
        bins3 = lk::createBinning(axis3);
    tos::AbsTwoPointGridPtr grid;
    if(polar) { 
        grid.reset(new tos::PolarGrid(bins1, bins2, bins3)); 
    } 
    else if (cart) { 
        std::cerr << "Not implemented yet!" << std::endl;
        return -1;
        // grid.reset(new tos::CartesianGrid(bins1, bins2, bins3)); 
    } 
    else {
        std::cerr << "Not implemented yet!" << std::endl;
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
    try {
        healxi(healbins, sight_lines, grid, max_ang, xi);
    }
    catch(std::exception const &e) {
        std::cerr << "Error while running the estimator: " << e.what() << std::endl;
    }

    // Save the estimator results
    try {
        std::vector<double> binCenters(3);
        std::ofstream out(outfile.c_str());
        for(int index = 0; index < xi.size(); ++index) {
            grid->getBinCenters(index, binCenters);
            out << index << ' ' << binCenters[0] << ' ' << binCenters[1] << ' ' << binCenters[2] << ' ' 
                << xi[index].didj << ' ' << xi[index].di << ' ' << xi[index].dj << ' ' << xi[index].wgt << std::endl;
        }
        out.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while saving results: " << e.what() << std::endl;
    }

    return 0;
}