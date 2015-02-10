// Created 5-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HDF5_DELTA
#define TURBOOCTOSPICE_HDF5_DELTA

#include "H5Cpp.h"
#include "types.h"

#include <vector>
#include <string>

namespace turbooctospice {

    class HDF5Delta {
    public:
        HDF5Delta(std::string filename);
        std::vector<Forest> loadForests(int ncombine, float forestlo, float foresthi, float speclo);
    private:
        std::string _filename;
    }; // HDF5Delta

} // turbooctospice

#endif