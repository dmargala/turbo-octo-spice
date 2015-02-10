// Created 24-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Read BOSS mock data 

#include "boost/program_options.hpp"
#include "boost/foreach.hpp"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include <iostream>
#include <fstream>

#include "turbooctospice.h"


namespace po = boost::program_options;
namespace lk = likely;
namespace tos = turbooctospice;

int main(int argc, char **argv) {

    std::string infile;
    po::options_description cli("Read mock data");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
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

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
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
        return -3;
    }
    if(verbose) {
        std::cout << "Read " << targetlist.size() << " rows from " << infile << std::endl;
    }

    long npixels(0);
    for(int i = 0; i < targetlist.size(); ++i){
        tos::MockSpectrum spectrum(targetlist[i], verbose);
        auto forest = spectrum.getForest();
        npixels += forest.pixels.size();
    }

    std::cout << "Number of pixels: " << npixels << std::endl;

    return 0;
}