// Created 5-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "HDF5Delta.h"
#include "H5Cpp.h"
#include "RuntimeError.h"

#include <cmath>
#include <algorithm>
#include <iostream>

namespace local = turbooctospice;

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

std::vector<local::Forest> local::HDF5Delta::loadForests() {
    // Open delta field file
    H5::Group grp;
    try {
        H5::H5File file(_filename, H5F_ACC_RDONLY);
        grp = file.openGroup("lines_of_sight");
    }
    catch(H5::FileIException const &e){
        throw local::RuntimeError("HDF5Delta: Could not open the specified file.");
        return -1;
    }

    static int nforests(0);
    struct kludge {
        static herr_t load_delta_los(hid_t loc_id, const char *name, const H5L_info_t *linfo, void *opdata) {
            H5::Group targetGroup;
            targetGroup = H5Gopen2(loc_id, name, H5P_DEFAULT);
            double ra, dec, z;
            readAttribute(targetGroup, "ra", ra);
            readAttribute(targetGroup, "dec", dec);
            readAttribute(targetGroup, "z", z);
            // Display progress
            if((nforests % 10000) == 0) std::cout << nforests << " : " << name << std::endl;
            // Read datasets
            std::vector<float> delta, r_comov, weight, loglam;
            readDataSet(targetGroup, "delta", delta);
            readDataSet(targetGroup, "weight", weight);
            readDataSet(targetGroup, "r_comov", r_comov);
            readDataSet(targetGroup, "loglam", loglam);
            // Init forest pixels
            Forest forest(ra, dec, nforests);
            // Iterate over pixels
            for(int i = 0; i < delta.size(); ++i) {
                // Save the pixel
                forest.pixels.push_back( {delta[i], loglam[i], weight[i], r_comov[i]} );
            }
            // Save forest
            ((std::vector<Forest>*) opdata)->push_back(forest);
            nforests++;
            targetGroup.close();
            return 0;
        }
    };

    // Load forests
    std::vector<Forest> forests;
    herr_t idx = H5Literate(grp.getId(), H5_INDEX_NAME, H5_ITER_INC, NULL, kludge::load_delta_los, &forests);
    return forests;
}
