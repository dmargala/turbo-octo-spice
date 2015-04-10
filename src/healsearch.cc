// Created 28-Jan-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "boost/program_options.hpp"
#include "boost/timer.hpp"

#include "turbooctospice.h"

namespace po = boost::program_options;
namespace tos = turbooctospice;

const double PI = std::atan(1.0)*4;
const double DEG2RAD = PI/180.0;

struct AngularPosition {
    double theta, phi, sth, cth, sph, cph;

    AngularPosition(double theta_, double phi_) : theta(theta_), phi(phi_) {
        sth = sin(theta);
        cth = cos(theta);
        sph = sin(phi);
        cph = cos(phi);
    }
};

double angularSeparation(double cth1, double sth1, double cph1, double sph1, double cth2, double sth2, double cph2, double sph2) {
    return sth1*sth2 + cth1*cth2*(sph1*sph2 + cph1*cph2);
}

double angularSeparation(AngularPosition const &q1, AngularPosition const &q2) {
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


void healxi(tos::HealpixBinsI const &healbins, std::vector<AngularPosition> const &forests, double maxAng) {

    double cosmin(std::cos(maxAng));
    std::vector<int> neighbors;
    long paircounts(0), closepairs(0);

    for(int i = 0; i < forests.size(); ++i) {
        auto qi = forests[i];
        // std::cout << i << " " << maxAng << std::endl;
        auto neighbors = healbins.getBinIndicesWithinRadius(qi.theta, qi.phi, maxAng);
        // for(int neighbor : neighbors) {
        //     std::cout << neighbor << " ";
        // }
        // std::cout << std::endl;
        // search neighboring healpix bins
        for(int neighbor : neighbors) {
            if(!healbins.checkBinExists(neighbor)) continue;
            // look for potential los pairs inside this healpix bin
            for(int j : healbins.getBinContents(neighbor)) {
                auto qj = forests[j];
                ++paircounts;
                // check angular separation
                double cosij(angularSeparation(qi, qj));
                // std::cout << "  " << j << " " << std::acos(cosij) << std::endl;
                if(cosij <= cosmin) continue;
                ++closepairs;
            }
        }
    }
    std::cout << "Number of pairs found: " << paircounts << std::endl;
    std::cout << "Number of close pairs found: " << closepairs << std::endl;
}

void dumbSearch(std::vector<AngularPosition> const &forests, double maxAng) {

    long paircounts(0), closepairs(0);
    double cosmin(std::cos(maxAng));

    for(auto qi : forests) {
        // std::cout << qi.theta << " " << qi.phi << " " << maxAng << std::endl;
        for(auto qj : forests) {
            ++paircounts;
            double cosij(angularSeparation(qi, qj));
            // std::cout << " " << qj.theta << " " << qj.phi << " " << maxAng << " " << std::acos(cosij) << std::endl;
            if(cosij <= cosmin) continue;
            ++closepairs;
        }
    }
    std::cout << "Number of pairs found: " << paircounts << std::endl;
    std::cout << "Number of close pairs found: " << closepairs << std::endl;
}

int main(int argc, char **argv) {

    int order, limit;
    double maxAng;
    std::string infile, outfile, axis1, axis2, axis3;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ("limit,n",po::value<int>(&limit)->default_value(1),
            "number of sightlines to search")
        ("order", po::value<int>(&order)->default_value(5),
            "Healpix map order parameter")
        ("radius", po::value<double>(&maxAng)->default_value(0.058),
            "Angular radius to search")
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

    boost::timer timer;

    // initialize Healpix bins
    tos::HealpixBinsI healbins(order);

    // load forests
    std::vector<AngularPosition> forests;// = {
    //     {0, 0},
    //     {0, .01},
    //     {0, .1},
    //     {0, .11},
    //     {0, .2},
    //     {0, .3},
    //     {0, .5}
    // };
    auto targetlist = readTargetList(infile);

    if(limit == 0 || limit > targetlist.size()) {
        limit = targetlist.size();
    }

    for(int i = 0; i < limit; ++i){
        tos::MockSpectrum spectrum(targetlist[i], verbose);
        forests.push_back(AngularPosition((90.0-spectrum.getDec())*DEG2RAD, spectrum.getRA()*DEG2RAD));
    }
    std::cout << "Finished loading forests: " << timer.elapsed() << std::endl;


    // add sight lines to healpix bins and calculate distances to pixels
    timer.restart();
    for(int i = 0; i < forests.size(); ++i) {
        healbins.addItem(forests[i].theta, forests[i].phi, i);
    }
    int numHealBinsOccupied(healbins.getNBins());
    std::cout << "Number of Healpix bins occupied: " << numHealBinsOccupied << " (" << static_cast<double>(numHealBinsOccupied)/(12*std::pow(4, order)) << ")" << std::endl;
    std::cout << "Finished filling healpix bins: " << timer.elapsed() << std::endl;

    std::cout << "Number of items added to healbins: " << healbins.getNEntries() << std::endl;
    // Generate the correlation function grid and run the estimator
    timer.restart();
    healxi(healbins, forests, maxAng);
    std::cout << "Finished finding pairs (healbins): " << timer.elapsed() << std::endl;

    timer.restart();
    dumbSearch(forests, maxAng);
    std::cout << "Finished finding pairs (dumb): " << timer.elapsed() << std::endl;

    return 0;
}