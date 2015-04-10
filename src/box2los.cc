// Created 23-Dec-2011 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>
// A tri-cubic 3D interpolation algorithm test program.

#include "H5Cpp.h"

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

void writeFloatAttr(H5::Group &group, std::string const &name, float value) {
    try {
        H5::Attribute attribute = group.createAttribute(name, H5::PredType::NATIVE_FLOAT, H5S_SCALAR);
        attribute.write(H5::PredType::NATIVE_FLOAT, &value);
    }
    catch(H5::AttributeIException const &e) {
        std::cerr << "Problem writing attribute : " << name << std::endl;
        e.printError();
    }
}

void writeFloatDSet(H5::Group &group, std::string const &name, std::vector<float> const &values) {
    try {
        const int rank = 1;
        hsize_t dims[rank] = { values.size() };
        H5::DataSet dataset = group.createDataSet(name, H5::PredType::NATIVE_FLOAT, {rank, dims});
        dataset.write(values.data(), H5::PredType::NATIVE_FLOAT);
    }
    catch(H5::DataSetIException const &e) {
        std::cerr << "Problem writing dataset : " << name << std::endl;
        e.printError();
    }
}

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

    H5::H5File file(outfile, H5F_ACC_TRUNC);

    H5::Group group(file.createGroup("/delta_field"));
    
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
        for(int trial = 0; trial < nlos; ++trial) {
            double phi = piovertwo*random->getUniform();
            double theta = std::acos(random->getUniform());
            double x = std::cos(phi)*std::sin(theta);
            double y = std::sin(phi)*std::sin(theta);
            double z = std::cos(theta);

            H5::Group target = group.createGroup(std::to_string(trial));

            writeFloatAttr(target, "ra", rad2deg*phi);
            writeFloatAttr(target, "dec", rad2deg*(piovertwo-theta));
            writeFloatAttr(target, "z", redshift);

            std::vector<float> delta, distance;
            for(int i = 0; i < nforest; ++i) {
                double r = rmin + i*dr;
                double d = interpolator(r*x, r*y, r*z);
                delta.push_back(d);
                distance.push_back(r);
                stats.accumulate(d);
            }

            writeFloatDSet(target, "delta", delta);
            writeFloatDSet(target, "distance", distance);

            target.close();
        }
        std::cout << "mean delta = " << stats.mean() << std::endl;
        std::cout << "RMS delta = " << stats.error() << std::endl;
    }
    catch(lk::RuntimeError const &e) {
        std::cerr << e.what() << std::endl;
    }
    catch(H5::GroupIException const &e) {
        e.printError();
    }

    group.close();
    file.close();

    return 0;
}
