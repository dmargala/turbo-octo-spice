// Created 5-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HDF5_DELTA
#define TURBOOCTOSPICE_HDF5_DELTA

#include "H5Cpp.h"
#include "types.h"
#include "turbooctospice.h"

#include <vector>
#include <string>

namespace turbooctospice {
    /// Represents an HDF5 delta field file.
    class HDF5Delta {
    public:
        /// Creates a HDF5 delta field file object
        /// @param filename The name of the file
        HDF5Delta(std::string filename);
        /// Loads forest sightlines from file
        /// @param keep_ngc specifies whether or not to include sight lines from the northern galactic cap
        /// @param keep_sgc specifies whether or not to include sight lines from the southern galactic cap
        std::vector<Forest> loadForests(bool keep_ngc=true, bool keep_sgc=true);
        // std::vector<Forest> loadForests(int ncombine, float forestlo, float foresthi, float speclo, bool debug=false);
    private:
        std::string _filename;
    }; // HDF5Delta

} // turbooctospice

#endif