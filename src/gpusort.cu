// nvcc -o gpusort src/gpusort.cu -std=c++11 -Iturbooctospice/ -lboost_program_options -lturbooctospice -Lbuild/ -lcfitsio -lCCfits

#include <thrust/sort.h> 
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/transform.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "boost/program_options.hpp"
#include "boost/timer.hpp"

namespace po = boost::program_options;

#include "turbooctospice.h"

namespace tos = turbooctospice;

const double PI = std::atan(1.0)*4;
const double DEG2RAD = PI/180.0;

struct AngularPosition {
    double theta, phi, sth, cth, sph, cph;

    AngularPosition() {};
    AngularPosition(double theta_, double phi_) : theta(theta_), phi(phi_) {
        sth = sin(theta);
        cth = cos(theta);
        sph = sin(phi);
        cph = cos(phi);
    }
};

__host__ __device__ double angularSeparation(double cth1, double sth1, double cph1, double sph1, double cth2, double sth2, double cph2, double sph2) {
    return sph1*sph2 + cph1*cph2*(sth1*sth2 + cth1*cth2);
};

__host__ __device__ double angularSeparation(AngularPosition const &p1, AngularPosition const &p2) {
    return angularSeparation(p1.cth, p1.sth, p1.cph, p1.sph, p2.cth, p2.sth, p2.cph, p2.sph);
};

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

struct rel_dist {
    const AngularPosition _q;

    rel_dist(AngularPosition q) : _q(q) { }

    __host__ __device__ 
    double operator()(AngularPosition const &p){
        return acos(angularSeparation(p, _q));
    }
};

template <typename T> struct 
is_greater_than {
    T threshold;
    is_greater_than(T thres) { threshold = thres; }

    __host__ __device__ bool operator()(const T& x) {
        return x > threshold;
    }
};

template <typename T> struct 
is_less_than {
    T threshold;
    is_less_than(T thres) { threshold = thres; }

    __host__ __device__ bool operator()(const T& x) {
        return x < threshold;
    }
};

struct within_radius {
    const AngularPosition _q;
    const double radius;

    within_radius(AngularPosition q, double r) : _q(q), radius(r) { }

    __host__ __device__ 
    bool operator()(AngularPosition const &p){
        return acos(angularSeparation(p, _q)) < radius;
    }
};

int main(int argc, char **argv) {

    std::string infile;
    int limit;
    double maxAng;
    po::options_description cli("Correlation function estimator");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ("limit,n",po::value<int>(&limit)->default_value(1),
            "number of sightlines to search")
        ("host", "Perform calculation on host too.")
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
    bool verbose(vm.count("verbose")), host(vm.count("host"));

    boost::timer timer;

    double cosmin(std::cos(maxAng));

    // load positions
    auto targetlist = readTargetList(infile);
    thrust::host_vector<AngularPosition> h_vec;
    for(int i = 0; i < targetlist.size(); ++i){
        tos::MockSpectrum spectrum(targetlist[i], verbose);
        // forests.push_back(spectrum.getForest(1, 1040.0, 1200.0, 3650.0));
        h_vec.push_back(AngularPosition((90.0-spectrum.getDec())*DEG2RAD, spectrum.getRA()*DEG2RAD));
    }
    std::cout << "Finished loading forests: " << timer.elapsed() << std::endl;

    if(limit == 0) {
        limit = targetlist.size();
    }

    // transfer data to the device
    timer.restart();
    thrust::device_vector<AngularPosition> d_vec = h_vec;
    thrust::device_vector<double> dist_d_vec(d_vec.size());
    std::cout << "Finished transfer to device: " << timer.elapsed() << std::endl;

    timer.restart();
    long paircounts = 0;
    for(int i = 0; i < limit; ++i) {
        // calculate separations
        //thrust::transform(d_vec.begin(), d_vec.end(), dist_d_vec.begin(), rel_dist(h_vec[i])); 
        within_radius pred(h_vec[i], maxAng);
        paircounts += thrust::count_if(d_vec.begin(), d_vec.end(), pred);
    }
    std::cout << "Number of pairs found: " << paircounts << std::endl;
    std::cout << "Finished transform/sort on device: " << timer.elapsed() << std::endl;


    if(host) {
        timer.restart();
        thrust::host_vector<double> dist_h_vec(h_vec.size());
        paircounts = 0;
        for(int i = 0; i < limit; ++i) {
            thrust::transform(h_vec.begin(), h_vec.end(), dist_h_vec.begin(), rel_dist(h_vec[i])); 
            paircounts += thrust::count_if(dist_h_vec.begin(), dist_h_vec.end(), is_less_than<double>(maxAng));
        }
        std::cout << "Number of pairs found: " << paircounts << std::endl;
        std::cout << "Finished transform/sort on host: " << timer.elapsed() << std::endl;
    }

    return 0;
}