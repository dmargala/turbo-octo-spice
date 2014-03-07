// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BIN_XYZ_PAIR
#define TURBOOCTOSPICE_BIN_XYZ_PAIR

#include "types.h"

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

namespace turbooctospice {

	template <typename T>
	class BinXYZPair {
	public:
        typedef T PixelType;
        BinXYZPair(lk::BinnedGrid const &grid, bool rmu, double x1min, double x1max, double x2min, double x2max): 
        _grid(grid), _rmu(rmu), _x1min(x1min), _x1max(x1max), _x2min(x2min), _x2max(x2max) {
        };
        ~BinXYZPair() {};
		void binPair(PixelType const &first, PixelType const &second, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
			// Calculate separation
            double dx = first.x - second.x;
            double dy = first.y - second.y;
            double dz = first.z - second.z;
            std::vector<double> _separation(2);
            if(_rmu) {
                _separation[0] = std::sqrt(dx*dx+dy*dy+dz*dz);
                _separation[1] = std::fabs(dz/_separation[0]);
            }
            else {
                _separation[0] = std::fabs(dz);
                _separation[1] = std::sqrt(dx*dx+dy*dy);
            }
            // Check that separation is within range of interest
            if(_separation[0] < _x1min || _separation[0] >= _x1max) return;
            if(_separation[1] < _x2min || _separation[1] >= _x2max) return;
            //std::cout << "(" << dx << "," << dy << "," << dz << ") -> (" << separation[0] << "," << separation[0] << ") -> ";

            // Bin pair
            try {
	            int index = _grid.getIndex(_separation);
	            double wgt = first.w*second.w;
	            double prod = wgt*first.d*second.d;
	            dsum[index] += prod;
	            wsum[index] += wgt;
	            nused++;
				//std::cout << index << " => " << prod << " (" << wgt << ")" << std::endl;
            }
            catch(lk::RuntimeError const &e) {
                std::cerr << "no xi bin found for i,j = " << first.i << ',' << second.i << std::endl;
            }
		};
        int getNBinsTotal() const { return _grid.getNBinsTotal(); };
	private:
        lk::BinnedGrid _grid;
        bool _rmu;
        double _x1min, _x1max, _x2min, _x2max;
	}; // BinXYZPair

    typedef BinXYZPair<Pixel> Bin;

} // turbooctospice

#endif // TURBOOCTOSPICE_BIN_XYZ_PAIR
