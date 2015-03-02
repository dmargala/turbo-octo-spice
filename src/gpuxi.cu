// Created 08-Aug-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Correlation function estimator using GPU.

// Compile using:
// nvcc -m64 -arch=sm_20 -lboost_program_options -llikely -lcosmo gpuxi.cu -o gpuxi

// Example usage:
// time -p ./gpuxi -i /Users/daniel/Cosmo/LyAlpha/cosmo/build/delta.dat --verbose

// nvcc -m64 -arch=sm_20 -c gpu/brute.cu -ccbin /usr/bin/clang
// g++ -I/usr/local/cuda/include -O3 src/gpuxi.cc -c
// g++ -O3 -llikely -L/usr/local/cuda/lib -lcuda -lcudart brute.o gpuxi.o -o gpuxi 


#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>

#include "../gpu/brute.h"

//#include "likely/likely.h"
// #include "boost/program_options.hpp"


// namespace po = boost::program_options;
// namespace lk = likely;


int main(int argc, char **argv) {

    // Configure command-line option processing
    std::string infile("sample-data/delta-64-10.dat"),outfile("output/gpu-xi-64-10.txt"),axis1("[0:200]*50"),axis2("[0:200]*50");
    long chunksize(1<<12);
    // po::options_description cli("Correlation function estimator");
    // cli.add_options()
    //     ("help,h", "Prints this info and exits.")
    //     ("verbose", "Prints additional information.")
    //     ("input,i", po::value<std::string>(&infile)->default_value(""),
    //         "Filename to read field samples from")
    //     ("output,o", po::value<std::string>(&outfile)->default_value("xi.dat"),
    //         "Filename to write correlation function to")
    //     ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
    //         "Axis-1 binning")
    //     ("axis2", po::value<std::string>(&axis2)->default_value("[0:200]*50"),
    //         "Axis-2 binning")
    //     ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")
    //     ("chunksize", po::value<long>(&chunksize)->default_value(4096),
    //         "Number of chunks to split the dataset into.")
    //     ;

    // //do the command line parsing now
    // po::variables_map vm;
    // try {
    //     po::store(po::parse_command_line(argc, argv, cli), vm);
    //     po::notify(vm);
    // }
    // catch(std::exception const &e) {
    //     std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
    //     return -1;
    // }
    // if(vm.count("help")) {
    //     std::cout << cli << std::endl;
    //     return 1;
    // }
    // bool verbose(vm.count("verbose")),rmu(vm.count("rmu"));
    bool verbose(true);

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > columns(5);

    std::ifstream in(infile.c_str());
    std::string line;
    while (std::getline(in, line)) {
        std::istringstream iss(line);
        float x, y, z, delta, wgt;
        if (!(iss >> x >> y >> z >> delta >> wgt)) { break; } // error
        columns[0].push_back(x);
        columns[1].push_back(y);
        columns[2].push_back(z);
        columns[3].push_back(delta);
        columns[4].push_back(wgt);
    }

    // try {
    //     std::ifstream in(infile.c_str());
    //     lk::readVectors(in, columns);
    //     in.close();
    // }
    // catch(std::exception const &e) {
    //     std::cerr << "Error while reading " << infile << ": " << e.what() << std::endl;
    //     return -3;
    // }
    if(verbose) {
        std::cout << "Read " << columns[0].size() << " rows from " << infile
            << std::endl;
    }

    // Generate the correlation function grid and run the estimator
    std::vector<double> xi;
    try {
        //lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2);
        // double x1min(bins1->getBinLowEdge(0)), x1max(bins1->getBinHighEdge(bins1->getNBins()-1));
        // double x2min(bins2->getBinLowEdge(0)), x2max(bins2->getBinHighEdge(bins2->getNBins()-1));
        // //lk::BinnedGrid grid(bins1,bins2);
        // int x1nbins = bins1->getNBins();

        double x1min(0), x1max(200);
        int x1nbins(50);

        bruteGPU(columns,x1min,x1max,x1nbins,xi,chunksize);
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
