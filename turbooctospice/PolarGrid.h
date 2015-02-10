// Created 09-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_POLAR_GRID
#define TURBOOCTOSPICE_POLAR_GRID

#include "AbsTwoPointGrid.h"

namespace turbooctospice {

    class PolarGrid : public AbsTwoPointGrid {
    public:
        PolarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3);
        virtual ~PolarGrid();
        double minAngularScale(double scale) const;
        double maxAngularScale(double scale) const;
        bool getSeparation(ForestPixel const &a, ForestPixel const &b,
        double const &cosij, double const &thetaij, std::vector<double> &separation) const {
            double distSq = a.distance*a.distance + b.distance*b.distance - 2*a.distance*b.distance*cosij;
            if(distSq >= x1maxSq || distSq < x1minSq) return false;
            separation[0] = std::sqrt(distSq);
            separation[1] = std::fabs(a.distance-b.distance)/separation[0];
            if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
            separation[2] = 0.5*(a.wavelength+b.wavelength) - logLyA;
            if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
            return true;
        }
    private:
        double x1minSq, x1maxSq;
    }; // PolarGrid

    inline double PolarGrid::minAngularScale(double scale) const { return xmin[0]*std::sqrt(1.0-xmax[1]*xmax[1])/scale; }
    inline double PolarGrid::maxAngularScale(double scale) const { return xmax[0]*std::sqrt(1.0-xmin[1]*xmin[1])/scale; }

} // turbooctospice

#endif // TURBOOCTOSPICE_POLAR_GRID
