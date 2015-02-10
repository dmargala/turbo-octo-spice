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
        AbsTwoPointGrid(lk::AbsBinningCPtr axis1, lk::AbsBinningCPtr axis2,
            lk::AbsBinningCPtr axis3);
        virtual ~AbsTwoPointGrid();
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

}

#endif
