// Created 5-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HDF5_DELTA
#define TURBOOCTOSPICE_HDF5_DELTA

#include "types.h"

#include <vector>
#include <string>

namespace turbooctospice {
    /// Represents an HDF5 delta field file.
    class HDF5Delta {
    public:
        /// Creates a HDF5 delta field file object
        /// @param filename The name of the file
        HDF5Delta(std::string filename);
        /// Loads forest sightlines from file. 
        /// Reads 'ra', 'dec', 'z', and 'plate' attributes from each sightline.
        /// Reads 'delta', 'weight', 'loglam', and 'r_comov' dataset from each sightline.
        /// @param groupname The name of the group containing forest sightlines
        std::vector<Forest> loadForests(std::string groupname);
    private:
        std::string _filename;
    }; // HDF5Delta

} // turbooctospice

#endif
