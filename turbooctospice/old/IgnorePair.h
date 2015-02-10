// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_DO_NOTHING
#define TURBOOCTOSPICE_DO_NOTHING

#include "types.h"

#include "likely/likely.h"

namespace lk = likely;

namespace turbooctospice {

	template <class T>
	class IgnorePair {
	public:
		typedef T PairType;
		void binPair(PairType const &pair, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
			nused++;
		};
		int getNBins() const { return 0; };
	private:
	}; // IgnorePair

    typedef IgnorePair<XYZPixelPair> Ignore;
    typedef IgnorePair<AngPixelPair> IgnoreAng;

} // turbooctospice

#endif // TURBOOCTOSPICE_DO_NOTHING
