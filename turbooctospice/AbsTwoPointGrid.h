// Created 6-Feb-2015 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_ABS_TWO_POINT_GRID
#define TURBOOCTOSPICE_ABS_TWO_POINT_GRID

#include "likely/likely.h"

#include "types.h"
#include "constants.h"

namespace lk = likely;

namespace turbooctospice {

    class AbsTwoPointGrid {
    public:
        AbsTwoPointGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3) :
            _grid(axis1, axis2, axis3) {
                for(int axis = 0; axis < 3; ++axis) {
                    auto bins = _grid.getAxisBinning(axis);
                    xmin.push_back(bins->getBinLowEdge(0));
                    xmax.push_back(bins->getBinHighEdge(bins->getNBins()-1));
                }
            };
        virtual ~AbsTwoPointGrid() {};
        // bool checkSeparation(std::vector<double> const &values) const;
        int getNBinsTotal() const;
        int getIndexNoCheck(std::vector<double> const &values,  std::vector<int> &binIndices) const;
        void getBinCenters(int index, std::vector<double> &binCenters) const;
        virtual double minAngularScale(double scale) const = 0;
        virtual double maxAngularScale(double scale) const = 0;
        virtual bool getSeparation(ForestPixel const &a, ForestPixel const &b, 
            double const &cosij, double const &thetaij, std::vector<double> &separation) const = 0;
    protected:
        std::vector<double> xmin, xmax;
        enum { LOS = 0, TRANSVERSE = 1, REDSHIFT = 2 };
    private:
        lk::BinnedGrid _grid;
    };

    inline int AbsTwoPointGrid::getNBinsTotal() const { return _grid.getNBinsTotal(); }
    inline int AbsTwoPointGrid::getIndexNoCheck(std::vector<double> const &values,  std::vector<int> &binIndices) const { 
        return _grid.getIndexNoCheck(values, binIndices); 
    }
    inline void AbsTwoPointGrid::getBinCenters(int index, std::vector<double> &binCenters) const {
        return _grid.getBinCenters(index, binCenters);
    }
    // inline bool AbsTwoPointGrid::checkSeparation(std::vector<double> const &values) const {
    //     for(int axis = 1; axis < 3; ++axis) 
    //         if(values[axis] < xmin[axis] || values[axis] >= xmax[axis]) return false;
    //     return true;
    // }

    class CartesianGrid : public AbsTwoPointGrid {
    public:
        CartesianGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3) :
            AbsTwoPointGrid(axis1, axis2, axis3) {
                x1minSq = xmin[0]*xmin[0];
                x1maxSq = xmax[0]*xmax[0];
            };
        virtual ~CartesianGrid() {};
        double minAngularScale(double scale) const;
        double maxAngularScale(double scale) const;
        bool getSeparation(ForestPixel const &a, ForestPixel const &b, 
        double const &cosij, double const &thetaij, std::vector<double> &separation) const {
            double distSq = a.distance*a.distance + b.distance*b.distance - 2*a.distance*b.distance*cosij;
            if(distSq >= x1maxSq || distSq < x1minSq) return false;
            separation[0] = std::fabs(a.distance-b.distance);
            separation[1] = std::sqrt(distSq - separation[0]*separation[0]);
            // separation[1] = thetaij*cosmology->getTransverseComovingScale(separation[2]);
            if(separation[1] < xmin[1] || separation[1] >= xmax[1]) return false;
            separation[2] = 0.5*(a.wavelength+b.wavelength) - logLyA;
            if(separation[2] < xmin[2] || separation[2] >= xmax[2]) return false;
            return true;
        }
    private:
        double x1minSq, x1maxSq;
    };

    inline double CartesianGrid::minAngularScale(double scale) const { return xmin[1]/scale; }
    inline double CartesianGrid::maxAngularScale(double scale) const { return xmax[1]/scale; }

    class PolarGrid : public AbsTwoPointGrid {
    public:
        PolarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3) :
            AbsTwoPointGrid(axis1, axis2, axis3) {
                x1minSq = xmin[0]*xmin[0];
                x1maxSq = xmax[0]*xmax[0];
            };
        virtual ~PolarGrid() {};
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
    };

    inline double PolarGrid::minAngularScale(double scale) const { return xmin[0]*std::sqrt(1.0-xmax[1]*xmax[1])/scale; }
    inline double PolarGrid::maxAngularScale(double scale) const { return xmax[0]*std::sqrt(1.0-xmin[1]*xmin[1])/scale; }

    class QuasarGrid : public AbsTwoPointGrid {
    public:
        QuasarGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2, lk::AbsBinningCPtr axis3) :
            AbsTwoPointGrid(axis1, axis2, axis3) {};
        virtual ~QuasarGrid() {};
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
    };

    inline double QuasarGrid::minAngularScale(double scale) const { return xmin[1]/rad2arcmin; }
    inline double QuasarGrid::maxAngularScale(double scale) const { return xmax[1]/rad2arcmin; }

}

#endif
