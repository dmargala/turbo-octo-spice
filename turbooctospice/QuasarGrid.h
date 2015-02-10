// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_QUASAR_GRID
#define TURBOOCTOSPICE_QUASAR_GRID

#include "AbsTwoPointGrid.h"

namespace turbooctospice {

    class QuasarGrid : public AbsTwoPointGrid {
    public:
        QuasarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3);
        virtual ~QuasarGrid();
        double minAngularScale(double scale) const;
        double maxAngularScale(double scale) const;
        bool getSeparation(ForestPixel const &a, ForestPixel const &b,
        double const &cosij, double const &thetaij, std::vector<double> &separation) const {
            separation[0] = std::fabs(a.wavelength-b.wavelength);
            if(separation[0] >= xmax[0] || separation[0] < xmin[0]) return false;
            separation[1] = thetaij*rad2arcmin;
            // we don't need to check thetaij, we've already done this for the line of sights
            // if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
            separation[2] = 0.5*(a.wavelength+b.wavelength) - logLyA;
            if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
            return true;
        }
    private:
    }; // QuasarGrid

    inline double QuasarGrid::minAngularScale(double scale) const { return xmin[1]/rad2arcmin; }
    inline double QuasarGrid::maxAngularScale(double scale) const { return xmax[1]/rad2arcmin; }

} // turbooctospice

#endif // TURBOOCTOSPICE_QUASAR_GRID
