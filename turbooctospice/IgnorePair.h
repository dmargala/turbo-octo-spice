// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_DO_NOTHING
#define TURBOOCTOSPICE_DO_NOTHING

#include "types.h"

#include "likely/likely.h"

namespace lk = likely;

namespace turbooctospice {
	template <typename PixelType>
	class DoNothing {
	public:
		void binPair(PixelType const &first, PixelType const &second, lk::BinnedGrid const &grid, bool rmu,
		double x1min, double x1max, double x2min, double x2max, std::vector<double> &dsum, std::vector<double> &wsum) const {};
	private:
	}; // DoNothing

} // turbooctospice

#endif // TURBOOCTOSPICE_DO_NOTHING
