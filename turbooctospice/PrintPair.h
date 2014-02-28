// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_PRINT_PAIR
#define TURBOOCTOSPICE_PRINT_PAIR

#include "types.h"

namespace turbooctospice {
	class PrintPair {
	public:
		PrintPair();
		virtual ~PrintPair();

		void binPair(PixelPair const &pair) const;
	private:
	}; // PrintPair
} // turbooctospice

#endif // TURBOOCTOSPICE_PRINT_PAIR
