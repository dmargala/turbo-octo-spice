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

double angularSeparation(double cth1, double sth1, double cph1, double sph1, double cth2, double sth2, double cph2, double sph2) {
    return sth1*sth2 + cth1*cth2*(sph1*sph2 + cph1*cph2);
}

double angularSeparation(tos::Forest const &q1, tos::Forest const &q2) {
    return angularSeparation(q1.cth, q1.sth, q1.cph, q1.sph, q2.cth, q2.sth, q2.cph, q2.sph);
}

std::vector<std::string> readTargetList(std::string const &infile) {
    // Read the input file
    std::vector<std::string> targetlist;
    try {
        std::ifstream in(infile.c_str());
        std::string line;
        while (std::getline(in, line)) {
            targetlist.push_back(line);
        }
        in.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while reading " << infile << ": " << e.what() << std::endl;
    }
    return targetlist;
}

struct XiEntry {
    double didj, di, dj, wgt;
    XiEntry() : didj(0), di(0), dj(0), wgt(0) {};
    XiEntry(double _di, double _dj, double _wgt) : didj(_di*_dj), di(_di), dj(_dj), wgt(_wgt) {};
};

void healxi(tos::HealpixBinsI const &healbins, std::vector<tos::Forest> const &forests,
tos::AbsTwoPointGridPtr const &grid, double minAng, double maxAng,
std::vector<XiEntry> &xi, std::vector<XiEntry> &xibin) {

    unsigned long numLOS(forests.size());

    std::cout << "Number of los: " << numLOS << std::endl;

    // stat counters
    unsigned long numHealpixBinsSearched(0), numLOSPairs(0), numLOSPairsUsed(0), 
        numPixels(0), numPixelPairs(0), numPixelPairsUsed(0);

    // create internal accumulation vectors
    int nbins = grid->getNBinsTotal();
    //std::vector<double> dsum(nbins,0), wsum(nbins,0), disum(nbins,0), djsum(nbins,0), wisum(nbins,0);
    std::vector<XiEntry> xisum(nbins, {});

    std::vector<double> dsum(nbins,0), wsum(nbins,0);

    // to avoid calling trig functions inside loop
    double cosmin(std::cos(maxAng)), cosmax(std::cos(minAng));

    // allocate temporary vectors before loop
    std::vector<int> neighbors;
    int binIndex;
    double distSq, pi_dist_sq, pi_projection_times_two;
    float weight, dprod;
    std::vector<int> binIndices(3);
    std::vector<double> separation(3);

    double min_dist = 0;
    double min_distsq = min_dist*min_dist;
    double max_distsq = 200*200;
    double bin_width = 4;

    boost::progress_display show_progress( forests.size() );

    for(int i = 0; i < forests.size(); ++i) {
        ++show_progress;
        // if((i % 10000) == 0) std::cout << i << std::endl;
        // if(xibin.size() > 2e7) break;
        auto qi = forests[i];
        numPixels += qi.pixels.size();
        //if(healbins.ang2pix(qi.p) != 4618) continue; // for debugging
        auto neighbors = healbins.getBinIndicesWithinRadius(qi.theta, qi.phi, maxAng);
        // search neighboring healpix bins
        for(int neighbor : neighbors) {
            ++numHealpixBinsSearched;
            if(!healbins.checkBinExists(neighbor)) continue;
            // look for potential los pairs inside this healpix bins
            for(int j : healbins.getBinContents(neighbor)) {
                // only count pairs once
                if(j <= i) continue;
                ++numLOSPairs;
                auto qj = forests[j];
                // check angular separation
                double cosij(angularSeparation(qi, qj));
                if(cosij <= cosmin || cosij >= cosmax) continue;
                double thetaij = std::acos(cosij);
                ++numLOSPairsUsed;
                // accumulate statistics for pixel pairs 
                for(auto &pi : qi.pixels) {
                    pi_dist_sq = pi.distance*pi.distance;
                    pi_projection_times_two = 2*pi.distance*cosij;
                    for(auto &pj : qj.pixels) {
                        ++numPixelPairs;
                        // check pairs are within our binning grid
                        distSq = pi_dist_sq + (pj.distance - pi_projection_times_two)*pj.distance;
                        if(distSq >= max_distsq || distSq < min_distsq) continue;
                        binIndex = int((std::sqrt(distSq) - min_dist)/bin_width);
                        //if(!grid->getBinIndex(pi, pj, cosij, thetaij, binIndex)) continue;
                        // if(!grid->checkSeparation(separation)) continue;
                        //
                        try {
                            weight = pi.weight*pj.weight;
                            dprod = weight*pi.value*pj.value;
                            // int binIndex = grid->getIndexNoCheck(separation, binIndices);
                            dsum[binIndex] += dprod;
                            // xisum[binIndex].di += pi.weight*pi.value;
                            // xisum[binIndex].dj += pj.weight*pj.value;
                            wsum[binIndex] += weight;
                            // if(index == nbins-1) {
                            //     xibin.push_back( XiEntry{pi.value, pj.value, 1.0} );
                            //     //xibin.push_back( XiEntry{pj.value, pi.value, 1.0} );
                            // }
                            ++numPixelPairsUsed;
                        }
                        catch(lk::RuntimeError const &e) {
                            std::cerr << "no xi bin found for i,j = " << i << ',' << j << std::endl;
                            std::cerr << separation[0] << " " <<  separation[1] << " " << separation[2] << std::endl;
                        }
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

    int order, combine;
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
        ("combine", po::value<int>(&combine)->default_value(10),
            "Number of wavelength bins to combine in fake spectra.")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:0.05]*25"),
            "Axis-1 binning, r_par (Mpc/h), r (Mpc/h), or dloglam")
        ("axis2", po::value<std::string>(&axis2)->default_value("[5:175]*1"),
            "Axis-2 binning, r_perp (Mpc/h), mu (r_par/r), or dtheta (arcmin)")
        ("axis3", po::value<std::string>(&axis3)->default_value("[0.4393:0.6284]*1"),
            "Axis-3 binning, redshift")
        ("polar", "(r,mu,z) binning")
        ("cart", "(r_perp,r_par,z) binning")
        ("fits", "read data from fits files, input must list of targets")
        ("debug", "use to specify debugging mode")
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
        debug(vm.count("debug"));

    // create grid for binning
    lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2), 
        bins3 = lk::createBinning(axis3);
    tos::AbsTwoPointGridPtr grid;
    if(polar) { 
        grid.reset(new tos::PolarGrid(bins1, bins2, bins3)); 
    } 
    else if (cart) { 
        grid.reset(new tos::CartesianGrid(bins1, bins2, bins3)); 
    } 
    else {
        grid.reset(new tos::QuasarGrid(bins1, bins2, bins3));
    }

    // set up cosmology
    cosmo::AbsHomogeneousUniversePtr cosmology;
    if(OmegaMatter == 0) OmegaMatter = 1 - OmegaLambda;
    cosmology.reset(new cosmo::LambdaCdmUniverse(OmegaLambda,OmegaMatter));

    // the minimum redshift sets the angular scale we will need to consider
    double zmin(std::pow(10, bins3->getBinLowEdge(0))-1);
    double scale(cosmology->getTransverseComovingScale(zmin));
    double minAng(grid->minAngularScale(scale)), maxAng(grid->maxAngularScale(scale));
    std::cout << "Transverse comoving scale at z = " << zmin <<  " (Mpc/h): " << scale << std::endl;
    std::cout << "Max (min) angular scale at z = " << zmin <<  " (rad): " << maxAng  << " (" << minAng << ")" << std::endl;

    // initialize Healpix bins
    tos::HealpixBinsI healbins(order);

    // load forests
    std::vector<tos::Forest> forests;
    if(fits) {
        auto targetlist = readTargetList(infile);
        for(int i = 0; i < targetlist.size(); ++i){
            tos::MockSpectrum spectrum(targetlist[i], verbose);
            forests.push_back(spectrum.getForest(combine, 1040.0, 1200.0, 3650.0));
        }
    }
    else {
        tos::HDF5Delta file(infile);
        forests = file.loadForests(combine, 1040.0, 1200.0, 3650.0, debug);
    }

    // add sight lines to healpix bins and calculate distances to pixels
    unsigned long totalpixels(0);
    for(int i = 0; i < forests.size(); ++i) {
        totalpixels += forests[i].pixels.size();
        healbins.addItem(forests[i].theta, forests[i].phi, i);
        if(!debug) {
            for(int j = 0; j < forests[i].pixels.size(); ++j) {
                float z(std::pow(10, forests[i].pixels[j].wavelength)/tos::lyA - 1.0);
                forests[i].pixels[j].distance = cosmology->getLineOfSightComovingDistance(z);
            }
        }
    }
    int numHealBinsOccupied(healbins.getNBins());
    std::cout << "Number of Healpix bins occupied: " << numHealBinsOccupied 
        << " (" << static_cast<double>(numHealBinsOccupied)/(12*std::pow(4, order)) << ")" << std::endl;
    std::cout << "Average number of pixels per LOS: " <<  static_cast<double>(totalpixels)/forests.size() << std::endl;

    // Generate the correlation function grid and run the estimator
    std::vector<XiEntry> xi, xibin;
    try {
        healxi(healbins, forests, grid, minAng, maxAng, xi, xibin);//, weights, meani, meanj);
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
            out << index << ' ' 
                << binCenters[0] << ' ' 
                << binCenters[1] << ' ' 
                << binCenters[2] << ' ' 
                << xi[index].didj << ' ' 
                << xi[index].di << ' '
                << xi[index].dj << ' '
                << xi[index].wgt << std::endl;
        }
        out.close();
    }
    catch(std::exception const &e) {
        std::cerr << "Error while saving results: " << e.what() << std::endl;
    }

    // Save the estimator results
    // try {
    //     std::ofstream out("xibin.txt");
    //     for(int index = 0; index < xibin.size(); ++index) {
    //         out << index << ' ' 
    //             << xibin[index].didj << ' ' 
    //             << xibin[index].di << ' '
    //             << xibin[index].dj << ' '
    //             << xibin[index].wgt << std::endl;
    //     }
    //     out.close();
    // }
    // catch(std::exception const &e) {
    //     std::cerr << "Error while saving results: " << e.what() << std::endl;
    // }

    return 0;
}