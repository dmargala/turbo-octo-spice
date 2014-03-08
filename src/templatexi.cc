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
    std::string infile,outfile,axis,buckets;

    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value(""),
            "Filename to write correlation function to")
        ("axis", po::value<std::string>(&axis)->default_value("[0:200]*50"),
            "Xi axis binning")
        ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")
        ("ignore", "No binnig (use for profiling pair search methods")
        ("buckets", po::value<std::string>(&buckets)->default_value(""),
            "bucket binning")
        ("norm", "Normalize xi by dividing by weights")
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
    bool verbose(vm.count("verbose")),rmu(vm.count("rmu")),ignore(vm.count("ignore")),norm(vm.count("norm"));

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

    // Parse the input file
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
    lk::AbsBinningCPtr bins = lk::createBinning(axis);
    int nbins(bins->getNBins());
    double min(bins->getBinLowEdge(0)), max(bins->getBinHighEdge(nbins-1));

    // Run the estimator
    std::vector<double> xi;

    if(buckets.length()) {
        lk::AbsBinningCPtr bucketbins1 = lk::createBinning(buckets),
            bucketbins2 = lk::createBinning(buckets),
            bucketbins3 = lk::createBinning(buckets);
        lk::BinnedGrid bucketgrid(bucketbins1, bucketbins2, bucketbins3);

        if(ignore) {
            tos::BucketIgnoreXi xiestimator(
                new tos::BucketSearch(pixels, bucketgrid, verbose), 
                new tos::Ignore, 
                verbose);
            xi = xiestimator.run(norm); 
        }
        else {
            tos::BucketBinXi xiestimator(
                new tos::BucketSearch(pixels, bucketgrid,  verbose), 
                new tos::Bin(min, max, nbins), 
                verbose);
            xi = xiestimator.run(norm); 
        }
    }
    else {
        if(ignore) {
            tos::BruteIgnoreXi xiestimator(
                new tos::BruteSearch(pixels, verbose), 
                new tos::Ignore, 
                verbose);
            xi = xiestimator.run(norm); 
        }
        else {
            tos::BruteBinXi xiestimator(
                new tos::BruteSearch(pixels, verbose), 
                new tos::Bin(min, max, nbins), 
                verbose);
            xi = xiestimator.run(norm); 
        }
    }

    // Save the estimator results
    if(outfile.length() > 0) {
        if(verbose) {
            std::cout << "Saving xi to " << outfile << std::endl;
        }
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