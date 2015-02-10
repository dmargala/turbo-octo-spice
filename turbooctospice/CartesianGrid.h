// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_CARTESIAN_GRID
#define TURBOOCTOSPICE_CARTESIAN_GRID

#include "AbsTwoPointGrid.h"

namespace turbooctospice {
    class CartesianGrid : public AbsTwoPointGrid {
    public:
        CartesianGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3);
        virtual ~CartesianGrid();
        double minAngularScale(double scale) const;
        double maxAngularScale(double scale) const;
        bool getSeparation(ForestPixel const &a, ForestPixel const &b, double const &cosij, double const &thetaij, std::vector<double> &separation) const;
    private:
        double x1minSq, x1maxSq;
	}; // CartesianGrid

    inline double CartesianGrid::minAngularScale(double scale) const { return xmin[1]/scale; }
    inline double CartesianGrid::maxAngularScale(double scale) const { return xmax[1]/scale; }

} // turbooctospice

#endif // TURBOOCTOSPICE_CARTESIAN_GRID
