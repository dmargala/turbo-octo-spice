// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_DO_NOTHING
#define TURBOOCTOSPICE_DO_NOTHING

#include "types.h"

#include "likely/likely.h"

namespace lk = likely;

namespace turbooctospice {
	template <typename PixelType>
	class IgnorePair {
	public:
		void binPair(PixelType const &first, PixelType const &second, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
			nused++;
		};
		int getNBinsTotal() const { return 0; };
	private:
	}; // IgnorePair

} // turbooctospice

#endif // TURBOOCTOSPICE_DO_NOTHING
