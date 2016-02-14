// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_POLAR_GRID
#define TURBOOCTOSPICE_POLAR_GRID

#include "AbsTwoPointGrid.h"

namespace turbooctospice {
    /// Represents a Polar grid for binning two point functions for LyA BAO analyses. 
    class PolarGrid : public AbsTwoPointGrid {
    public:
        enum AxisLabels {
            R, Mu, LogLya1pz
        };
        /// Creates a new grid with three axes.
        /// @param axis1 Line of sight separation axis binning : |r|
        /// @param axis2 Transverse separation axis binning : mu
        /// @param axis3 Mean distance axis binning : log(1+z)
        PolarGrid(likely::AbsBinningCPtr axis1, likely::AbsBinningCPtr axis2, likely::AbsBinningCPtr axis3);
        virtual ~PolarGrid();
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

        bool getBinIndex(ForestPixel const &a, ForestPixel const &b, 
            double const &cosij, double const &thetaij, int &binIndex) const;
    private:
        double x1minSq, x1maxSq;
    }; // PolarGrid

    inline double PolarGrid::minAngularScale(double scale) const { return xmin[R]*std::sqrt(1.0-xmax[Mu]*xmax[Mu])/scale; }
    inline double PolarGrid::maxAngularScale(double scale) const { return xmax[R]*std::sqrt(1.0-xmin[Mu]*xmin[Mu])/scale; }

} // turbooctospice

#endif // TURBOOCTOSPICE_POLAR_GRID
