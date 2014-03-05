// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_PRINT_PAIR
#define TURBOOCTOSPICE_PRINT_PAIR

#include "types.h"

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

namespace turbooctospice {
	template <typename PixelType>
	class XYZPair {
	public:
		void binPair(PixelType const &first, PixelType const &second, lk::BinnedGrid const &grid, bool rmu,
		double x1min, double x1max, double x2min, double x2max, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
			// Calculate separation
            double dx = first.x - second.x;
            double dy = first.y - second.y;
            double dz = first.z - second.z;
            std::vector<double> separation(2);
            if(rmu) {
                separation[0] = std::sqrt(dx*dx+dy*dy+dz*dz);
                separation[1] = std::fabs(dz/separation[0]);
            }
            else {
                separation[0] = std::fabs(dz);
                separation[1] = std::sqrt(dx*dx+dy*dy);
            }
            // Check that separation is within range of interest
            if(separation[0] < x1min || separation[0] >= x1max) return;
            if(separation[1] < x2min || separation[1] >= x2max) return;
            //std::cout << "(" << dx << "," << dy << "," << dz << ") -> (" << separation[0] << "," << separation[0] << ") -> ";

            // Bin pair
            try {
	            int index = grid.getIndex(separation);
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
	private:
	}; // XYZPair
} // turbooctospice

#endif // TURBOOCTOSPICE_PRINT_PAIR
