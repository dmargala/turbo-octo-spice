// Created 30-Jan-2013 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// KDTree test program for correlation function estimation

// g++ -02 -Wall -g -lkdtree -lboost_program_options -llikely -lcosmo -o kdtreetest kdtreetest.c
// time -p ./kdtreetest -i sample-data/delta-tinier.dat --verbose



#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include "kdtree.h"

#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/format.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>


namespace po = boost::program_options;
namespace lk = likely;


int main(int argc, char **argv) {

	// Configure command-line option processing
    std::string infile,outfile,axis1,axis2;
    float rmax;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value("xi.dat"),
            "Filename to write correlation function to")
        ("rmax", po::value<float>(&rmax)->default_value(200.),
        	"Maximum distance for pair search")
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

	kdtree *kd;
	kd = kd_create(3);

	long n = columns[0].size();

	for(long i=0; i < n; ++i) {
		float x, y, z;
		x = columns[0][i];
		y = columns[1][i];
		z = columns[2][i];

		assert(kd_insert3(kd, x, y, z, 0) == 0);
	}
	if(verbose) {
		std::cout << "Finished building tree." << std::endl;
	}

	unsigned long npairs = 0;

	for(long i=0; i < n; ++i) {
		if(verbose) {
			if(i % 1000 == 0){
				std::cout << "Working on " << i << " out of " << n << " points." << std::endl;
			}
		}
		float x, y, z;
		x = columns[0][i];
		y = columns[1][i];
		z = columns[2][i];

		kdres *set;
		set = kd_nearest_range3(kd, x, y, z, rmax);
		npairs += kd_res_size(set);
		kd_res_free(set);
	}

	std::cout << "Range query returned " << npairs << " pairs out of " << n*(n-1)/2 << " possible." << std::endl;

	kd_free(kd);
	return 0;
}
