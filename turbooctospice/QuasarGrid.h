// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_QUASAR_GRID
#define TURBOOCTOSPICE_QUASAR_GRID

#include "AbsTwoPointGrid.h"

namespace turbooctospice {
    /// Represents a Quasar grid for binning two point functions for LyA BAO analyses. 
    class QuasarGrid : public AbsTwoPointGrid {
    public:
        /// Creates a new grid with three axes.
        /// @param axis1 Line of sight separation axis binning : log(lam1/lam2)
        /// @param axis2 Transverse separation axis binning : angular separation (arcminutes)
        /// @param axis3 Mean distance axis binning : log(1+z)
        QuasarGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3);
        virtual ~QuasarGrid();
        /// Returns the minimum angular scale of the grid in radians from the provided
        /// transervse comving scale.
        /// @param scale Transverse comoving distance (Mpc/h) at minimum redshift
        double minAngularScale(double scale) const;
        /// Returns the maximum angular scale of the grid in radians from the provided
        /// transervse comving scale.
        /// @param scale Transverse comoving distance (Mpc/h) at minimum redshift
        double maxAngularScale(double scale) const;
        /// Tests if the pair of pixels are binable and fills the vector provided 
        /// with coordinate values along each axis.
        /// @param a First ForestPixel of pair
        /// @param b Second ForestPixel of pair
        /// @param cosij Cosine of the angular separation between the line of sights to each pixel
        /// @param thetaij Angular separation between the line of sights to each pixel
        /// @param separation Output vector containing the separation of the pixel pairs along each axis
        bool getSeparation(ForestPixel const &a, ForestPixel const &b,
            double const &cosij, double const &thetaij, std::vector<double> &separation) const;
    private:
    }; // QuasarGrid

    inline double QuasarGrid::minAngularScale(double scale) const { return xmin[1]/rad2arcmin; }
    inline double QuasarGrid::maxAngularScale(double scale) const { return xmax[1]/rad2arcmin; }

} // turbooctospice

#endif // TURBOOCTOSPICE_QUASAR_GRID
