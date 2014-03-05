// Created 27-Feb-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Example strategy (polcicy) design pattern using templated classes

#include <fstream>

#include "boost/program_options.hpp"

#include "turbooctospice.h"
#include "likely/likely.h"

namespace lk = likely;

namespace tos = turbooctospice;
namespace po = boost::program_options;

int main(int argc, char **argv) {

    // Configure command-line option processing
    std::string infile,outfile,axis1,axis2;

    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value(""),
            "Filename to write correlation function to")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
            "Axis-1 binning")
        ("axis2", po::value<std::string>(&axis2)->default_value("[0:200]*50"),
            "Axis-2 binning")
        ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")
        ("ignore", "No binnig (use for profiling pair search methods")
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
    bool verbose(vm.count("verbose")),rmu(vm.count("rmu")),ignore(vm.count("ignore"));

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }
    std::vector<std::vector<double> > columns(5);
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

    tos::Pixels pixels;
    for(int i = 0; i < columns[0].size(); ++i) {
        tos::Pixel pixel;
        pixel.x = columns[0][i];
        pixel.y = columns[1][i];
        pixel.z = columns[2][i];
        pixel.d = columns[3][i];
        pixel.w = columns[4][i];
        pixel.i = i;
        pixels.push_back(pixel);
    }

    // Instantiate the correlation function grid
    lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2);
        double x1min(bins1->getBinLowEdge(0)), x1max(bins1->getBinHighEdge(bins1->getNBins()-1));
        double x2min(bins2->getBinLowEdge(0)), x2max(bins2->getBinHighEdge(bins2->getNBins()-1));
    lk::BinnedGrid grid(bins1,bins2);

    // Run the estimator
    std::vector<double> xi;
 
    typedef tos::BrutePairSearch<tos::Pixels> BPS;
    typedef tos::IgnorePair<tos::Pixel> IP;
    typedef tos::BinXYZPair<tos::Pixel> BP;
    typedef tos::XiEstimator<BPS, IP> IgnoreXi;
    typedef tos::XiEstimator<BPS, BP> BruteXi;

    if(ignore) {
        IgnoreXi ignorexi(new BPS, new IP);
        ignorexi.run(pixels, pixels, xi); 
    }
    else {
        BruteXi brutexi(new BPS, new BP(grid, rmu, x1min, x1max, x2min, x2max));
        brutexi.run(pixels, pixels, xi); 
    }

    long n = columns[0].size();
    std::cout << "Number of distinct pairs " << n*(n-1)/2 << std::endl;

    // Save the estimator results
    if(0 == outfile.length()) {
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