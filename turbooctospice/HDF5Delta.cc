// Created 5-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "HDF5Delta.h"
#include "H5Cpp.h"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace local = turbooctospice;

const double PI = std::atan(1.0)*4;
const double DEG2RAD = PI/180.0;
const double lyA = 1216.0;

template<typename T> void readDataSet(H5::Group& grp, const std::string& name, std::vector<T>& data) {
    // Open the data set
    H5::DataSet dataset = grp.openDataSet(name);
    // Need to figure out dataset dimensions
    H5::DataSpace dataspace = dataset.getSpace();
    // Get the number of dimensions in the dataspace.
    // int rank = dataspace.getSimpleExtentNdims();
    // std::cout << "\tndims: " << rank << std::endl;
    hsize_t dims_out[1];
    dataspace.getSimpleExtentDims(dims_out, NULL);
    int size = dims_out[0];
    // create a vector the same size as the dataset
    data.resize(size);
    // read dataset into vector
    dataset.read(data.data(), H5::PredType::NATIVE_FLOAT);
    // close dataset
    dataset.close();
}

template<typename T> void readAttribute(H5::Group& grp, const std::string& name, T& value){
    H5::Attribute attr;
    attr = grp.openAttribute(name);
    attr.read(attr.getDataType(), &value);
}

local::HDF5Delta::HDF5Delta(std::string filename) : _filename(filename) {};

// Operator function
//extern "C" herr_t load_quasar(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata);

std::vector<local::Forest> local::HDF5Delta::loadForests(int ncombine, float forestlo, float foresthi, float speclo, bool debug) {
    // open delta field file
    H5::H5File file(_filename, H5F_ACC_RDONLY);
    H5::Group grp = file.openGroup("lines_of_sight");

    static float flo(forestlo), fhi(foresthi), slo(speclo);
    static int n(ncombine), nforests(0);
    static bool select_south_cap(debug);

    struct kludge { 
        static herr_t load_quasar(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
            H5::Group targetGroup;
            targetGroup = H5Gopen2(loc_id, name, H5P_DEFAULT);

            if((nforests % 10000) == 0) std::cout << nforests << " : " << name << std::endl;

            double ra, dec, z;
            readAttribute(targetGroup, "ra", ra);
            readAttribute(targetGroup, "dec", dec);
            readAttribute(targetGroup, "z", z);

            // Read datasets
            std::vector<float> redshifts, delta, ivar;
            readDataSet(targetGroup, "absorber_z", redshifts);
            readDataSet(targetGroup, "absorber_delta", delta);
            readDataSet(targetGroup, "absorber_ivar", ivar);

            // init forest pixels
            Forest forest;
            double theta((90.0-dec)*DEG2RAD), phi(ra*DEG2RAD);
            forest.phi = phi;
            forest.theta = theta;
            forest.sth = std::sin(theta);
            forest.cth = std::cos(theta);
            forest.sph = std::sin(phi);
            forest.cph = std::cos(phi);

            // Max and min wavelength defined by lyman alpha forest range and spec cutoff
            float minLogLambda(std::log10(((flo*(1+z) > slo) ? flo*(1+z) : slo)));
            float maxLogLambda(std::log10(fhi*(1+z)));
            // Iterate over pixels
            float zavg, davg, wavg, distance, logLambda;
            for(int i = 0; i < redshifts.size()-n+1; i += n) {
                zavg = 0; davg = 0; wavg = 0;
                for(int j = 0; j < n; ++j) {
                    zavg += redshifts[i+j];
                    davg += delta[i+j];
                    wavg += ivar[i+j];
                    distance += 0; // ???
                }
                zavg /= n;
                logLambda = std::log10(lyA*(1+zavg));
                // If below forest, skip to next pixel
                if (logLambda < minLogLambda) continue;
                // If beyond forest, we're done
                if (logLambda > maxLogLambda) break;
                // Save the pixel
                forest.pixels.push_back( {davg/n, logLambda, wavg/n, distance/n} );
            }

            ((std::vector<Forest>*) opdata)->push_back(forest);
            nforests++;

            targetGroup.close();

            return 0;
        }

        static herr_t load_delta_los(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
            H5::Group targetGroup;
            targetGroup = H5Gopen2(loc_id, name, H5P_DEFAULT);

            double ra, dec, z;
            readAttribute(targetGroup, "ra", ra);
            readAttribute(targetGroup, "dec", dec);
            readAttribute(targetGroup, "z", z);

            ra *= 180.0/PI;
            dec *= 180.0/PI;

            if(select_south_cap && !(ra < 120 || ra > 330)) {
                if((nforests % 10000) == 0) std::cout << nforests << " : " << name << std::endl;
    
                // Read datasets
                std::vector<float> delta, r_comov, weight, loglam;
                readDataSet(targetGroup, "delta", delta);
                readDataSet(targetGroup, "weight", weight);
                readDataSet(targetGroup, "r_comov", r_comov);
                readDataSet(targetGroup, "loglam", loglam);

                // init forest pixels
                Forest forest;
                double theta((90.0-dec)*DEG2RAD), phi(ra*DEG2RAD);
                forest.phi = phi;
                forest.theta = theta;
                forest.sth = std::sin(theta);
                forest.cth = std::cos(theta);
                forest.sph = std::sin(phi);
                forest.cph = std::cos(phi);

                // Iterate over pixels
                for(int i = 0; i < delta.size(); ++i) {
                    // Save the pixel
                    forest.pixels.push_back( {delta[i], loglam[i], weight[i], r_comov[i]} );
                }

                ((std::vector<Forest>*) opdata)->push_back(forest);
                nforests++;
            }

            targetGroup.close();

            return 0;
        }

    };

    std::vector<Forest> forests;

    // Load data
    if(debug) {
        herr_t idx = H5Literate(grp.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, kludge::load_delta_los, &forests);
    }
    else {
        herr_t idx = H5Literate(grp.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, kludge::load_quasar, &forests);
    }

    return forests;
}

