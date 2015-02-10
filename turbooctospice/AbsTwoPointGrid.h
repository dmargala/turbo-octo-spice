// Created 6-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_ABS_TWO_POINT_GRID
#define TURBOOCTOSPICE_ABS_TWO_POINT_GRID

#include "likely/likely.h"

#include "types.h"
#include "constants.h"

namespace turbooctospice {
    /// Represents a grid for binning two point functions for LyA BAO analyses. 
    class AbsTwoPointGrid {
    public:
        /// Creates a new grid with three axes.
        /// @param axis1 Line of sight separation axis binning
        /// @param axis2 Transverse separation axis binning
        /// @param axis3 Mean distance axis binning
        AbsTwoPointGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2,
            likely::AbsBinningCPtr axis3);
        virtual ~AbsTwoPointGrid();
        /// Returns the total number of bins covering the rectangular volume of this grid.
        int getNBinsTotal() const;
        /// Returns the global index corresponding to the specified coordinate values along
        /// each axis.
        /// @param values Coordinate values.
        /// @param binIndices Temporary vector used to store bin indices.
        int getIndexNoCheck(std::vector<double> const &values,  std::vector<int> &binIndices) const;
        /// Fills the vector provided with the bin centers along each axis for the specified
        /// global index.
        /// @param index Global bin index.
        /// @param binCenters Output vector containing bin centers for the specified bin index.
        void getBinCenters(int index, std::vector<double> &binCenters) const;
        /// Returns the minimum angular scale of the grid in radians from the provided
        /// transervse comving scale.
        /// @param scale Transverse comoving distance (Mpc/h) at minimum redshift
        virtual double minAngularScale(double scale) const = 0;
        /// Returns the maximum angular scale of the grid in radians from the provided
        /// transervse comving scale.
        /// @param scale Transverse comoving distance (Mpc/h) at minimum redshift
        virtual double maxAngularScale(double scale) const = 0;
        /// Tests if the pair of pixels are binable and fills the vector provided 
        /// with coordinate values along each axis.
        /// @param a First ForestPixel of pair
        /// @param b Second ForestPixel of pair
        /// @param cosij Cosine of the angular separation between the line of sights to each pixel
        /// @param thetaij Angular separation between the line of sights to each pixel
        /// @param separation Output vector containing the separation of the pixel pairs along each axis
        virtual bool getSeparation(ForestPixel const &a, ForestPixel const &b, 
            double const &cosij, double const &thetaij, std::vector<double> &separation) const = 0;
    protected:
        std::vector<double> xmin, xmax;
        //enum { LOS = 0, TRANSVERSE = 1, REDSHIFT = 2 };
    private:
        likely::BinnedGrid _grid;
    };

    inline int AbsTwoPointGrid::getNBinsTotal() const { return _grid.getNBinsTotal(); }
    inline int AbsTwoPointGrid::getIndexNoCheck(std::vector<double> const &values,  std::vector<int> &binIndices) const { 
        return _grid.getIndexNoCheck(values, binIndices); 
    }
    inline void AbsTwoPointGrid::getBinCenters(int index, std::vector<double> &binCenters) const {
        return _grid.getBinCenters(index, binCenters);
    }

}

#endif
