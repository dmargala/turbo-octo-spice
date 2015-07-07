// Created 23-Dec-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// A tri-cubic 3D interpolation algorithm test program.

#include "likely/likely.h"

#include "boost/program_options.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"

#include <iostream>
#include <fstream>
#include <cmath>

namespace lk = likely;
namespace po = boost::program_options;

const double piovertwo = std::atan2(1,0);
const double rad2deg = 90.0/piovertwo;

template <typename T>
struct FieldSample {
    T x, y, z, d, w;
    FieldSample(T _x, T _y, T _z, T _d, T _w) 
    : x(_x), y(_y), z(_z), d(_d), w(_w) {};
};

int main(int argc, char **argv) {

    // Configure command-line option processing
    int nx, ny, nz, nlos, nforest;
    double spacing, rmin, rmax, redshift;
    std::string infile, outfile;
    po::options_description cli("Tri-cubic interpolation test program");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read field samples from")
        ("output,o", po::value<std::string>(&outfile)->default_value(""),
            "Filename to write field line of sight samples to")
        ("nx", po::value<int>(&nx)->default_value(256),
            "Number of subdivisions along grid x axis.")
        ("ny", po::value<int>(&ny)->default_value(0),
            "Number of subdivisions along grid y axis.")
        ("nz", po::value<int>(&nz)->default_value(0),
            "Number of subdivisions along grid z axis.")
        ("spacing", po::value<double>(&spacing)->default_value(2),
            "Spacing between grid points")
        ("nlos", po::value<int>(&nlos)->default_value(1000),
            "Number of random points for testing the interpolator.")
        ("rmin", po::value<double>(&rmin)->default_value(4000),
            "Min forest distance")
        ("rmax", po::value<double>(&rmax)->default_value(4300),
            "Max forest distance")
        ("nforest", po::value<int>(&nforest)->default_value(100),
            "Number of forest pixels")
        ("redshift", po::value<double>(&redshift)->default_value(3),
            "Default quasar redshift")
        ("text","write to text file")
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

    // Fill in missing dimensions
    if(ny == 0) ny = nx;
    if(nz == 0) nz = nx;

    if(nx <= 0 || ny <= 0 || nz <= 0) {
        std::cerr << "Bad dimensions nx,ny,nz." << std::endl;
        return -1;
    }

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -1;
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

    try {
        int n(nx*ny*nz);
        lk::TriCubicInterpolator::DataCube data(new double[n]);
        for(int ix = 0; ix < nx; ++ix) {
            for(int iy = 0; iy < ny; ++iy) {
                for(int iz = 0; iz < nz; ++iz) {
                    data[ix + nx*(iy + ny*iz)] = columns[3][ix + nx*(iy + ny*iz)];
                }
            }
        }
        // Interpolate in this datacube at random points.
        lk::TriCubicInterpolator interpolator(data, spacing, nx, ny, nz);
        lk::RandomPtr random = lk::Random::instance();
        random->setSeed(1234);
        lk::WeightedAccumulator stats;
        double dr = (rmax-rmin)/(nforest+1);
        std::cout << "rmin, rmax, dr: " << rmin << ", " << rmax << ", " << dr << std::endl;
        std::vector<FieldSample<double> > samples;
        for(int trial = 0; trial < nlos; ++trial) {
            double phi = piovertwo*random->getUniform();
            double theta = std::acos(random->getUniform());
            double x = std::cos(phi)*std::sin(theta);
            double y = std::sin(phi)*std::sin(theta);
            double z = std::cos(theta);

            for(int i = 0; i < nforest; ++i) {
                double r = rmin + i*dr;
                double d = interpolator(r*x, r*y, r*z);
                samples.push_back({r*x, r*y, r*z, d, 1});
                stats.accumulate(d);
            }
        }
        std::cout << "mean delta = " << stats.mean() << std::endl;
        std::cout << "RMS delta = " << stats.error() << std::endl;

        // Save text file
        try {
            std::ofstream out(outfile.c_str());
            for(int i = 0; i < samples.size(); ++i) {
                out << samples[i].x << ' ' 
                    << samples[i].y << ' ' 
                    << samples[i].z << ' ' 
                    << samples[i].d << ' ' 
                    << samples[i].w << std::endl;
            }
            out.close();
        }
        catch(std::exception const &e) {
            std::cerr << "Error while saving results: " << e.what() << std::endl;
        }
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }


    return 0;
}
