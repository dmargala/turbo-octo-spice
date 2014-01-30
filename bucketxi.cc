// Created 30-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Bucket correlation function estimator.

// g++ -02 -lboost_program_options -llikely -lcosmo -o bucketxi bucketxi.cc 

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>

namespace po = boost::program_options;
namespace lk = likely;

void bucketXi(std::vector<std::vector<double> > const &columns, lk::BinnedGrid const &grid, bool rmu,
double x1min, double x1max, double x2min, double x2max, std::vector<double> &xi, lk::BinnedGrid const &bucketgrid) {
    std::cout << "Starting bucketxi..." << std::endl;
    // create internal accumulation vectors
    int nbins = grid.getNBinsTotal();
    std::vector<double> dsum(nbins,0.), wsum(nbins,0.);
    // Do a brute force loop over all pairs
    long n(columns[0].size());
    long npair(0), nused(0);

    std::map<int,std::vector<int> > buckets;

    std::vector<double> position(3);
    for(int i = 0; i < n; ++i) {
        position[0] = columns[0][i];
        position[1] = columns[1][i];
        position[2] = columns[2][i];
        int bucketIndex = bucketgrid.getIndex(position);
        if(buckets.count(bucketIndex) > 0) {
            buckets[bucketIndex].push_back(i);
        }
        else {
            buckets[bucketIndex] = std::vector<int>(1,i);
        }
    }

    int nbuckets = buckets.size();
    std::cout << "We have " << nbuckets << " buckets" << std::endl;
    std::vector<double> separation(2);

    std::map<int,std::vector<int> >::iterator iter;
    for(iter = buckets.begin(); iter != buckets.end(); ++iter) {
        std::vector<int> bucket = iter->second;  
        int npoints = bucket.size();
        for(int iBucket = 0; iBucket < npoints-1; ++iBucket) {
            int i = bucket[iBucket];
            double xi = columns[0][i];
            double yi = columns[1][i];
            double zi = columns[2][i];
            double di = columns[3][i];
            double wi = columns[4][i];
            for(int jBucket = iBucket+1; jBucket < npoints; ++jBucket) {
                int j = bucket[jBucket];
                double dx = xi - columns[0][j];
                double dy = yi - columns[1][j];
                double dz = zi - columns[2][j];
                if(rmu) {
                    separation[0] = std::sqrt(dx*dx+dy*dy+dz*dz);
                    separation[1] = std::fabs(dz/separation[0]);
                }
                else {
                    separation[0] = std::fabs(dz);
                    separation[1] = std::sqrt(dx*dx+dy*dy);
                }
                npair++;
                if(separation[0] < x1min || separation[0] >= x1max) continue;
                if(separation[1] < x2min || separation[1] >= x2max) continue;
                try {
                    int index = grid.getIndex(separation);
                    double wgt = wi*columns[4][j];
                    dsum[index] += wgt*di*columns[3][j];
                    wsum[index] += wgt;
                    nused++;
                }
                catch(lk::RuntimeError const &e) {
                    std::cerr << "no bin found for i,j = " << i << ',' << j << std::endl;
                    // std::cerr << ".. separation = (" << separation[0] <<  "," << separation[1] << ")" << std::endl;
                    // std::cerr << ".. (dx,dy,dz) = (" << dx << "," << dy << "," << dz << ")" << std::endl;
                }
            }
        }
    }

    std::cout << "used " << nused << " of " << npair << " pairs." << std::endl;
    std::cout << "Number of distinct pairs " << n*(n-1)/2 << std::endl;
    for(int index = 0; index < nbins; ++index) {
        if(wsum[index] > 0) dsum[index] /= wsum[index];
    }
    dsum.swap(xi);
}

int main(int argc, char **argv) {
    
    // Configure command-line option processing
    std::string infile,outfile,axis1,axis2,buckets;

    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value("xi.dat"),
            "Filename to write correlation function to")
        ("axis1", po::value<std::string>(&axis1)->default_value("[0:200]*50"),
            "Axis-1 binning")
        ("axis2", po::value<std::string>(&axis2)->default_value("[0:200]*50"),
            "Axis-2 binning")
        ("rmu", "Use (r,mu) binning instead of (rP,rT) binning")
        ("buckets", po::value<std::string>(&buckets)->default_value("[0:1000]*5"),
            "bucket binning")
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
    bool verbose(vm.count("verbose")),rmu(vm.count("rmu"));

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

    // Generate the correlation function grid and run the estimator
    std::vector<double> xi;
    try {
        lk::AbsBinningCPtr bins1 = lk::createBinning(axis1), bins2 = lk::createBinning(axis2);
        double x1min(bins1->getBinLowEdge(0)), x1max(bins1->getBinHighEdge(bins1->getNBins()-1));
        double x2min(bins2->getBinLowEdge(0)), x2max(bins2->getBinHighEdge(bins2->getNBins()-1));
        lk::BinnedGrid grid(bins1,bins2);

        lk::AbsBinningCPtr bucketbins1 = lk::createBinning(buckets),
            bucketbins2 = lk::createBinning(buckets),
            bucketbins3 = lk::createBinning(buckets);
        lk::BinnedGrid bucketgrid(bucketbins1, bucketbins2, bucketbins3);
        bucketXi(columns,grid,rmu,x1min,x1max,x2min,x2max,xi,bucketgrid);
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
