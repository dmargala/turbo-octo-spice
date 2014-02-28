// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_DO_NOTHING
#define TURBOOCTOSPICE_DO_NOTHING

#include "types.h"

namespace turbooctospice {
	class DoNothing {
	public:
		DoNothing() {};
		virtual ~DoNothing() {};

		void binPair(PixelPair const &pair) const {};
	private:
	}; // DoNothing

} // turbooctospice

#endif // TURBOOCTOSPICE_DO_NOTHING
