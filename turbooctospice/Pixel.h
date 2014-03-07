// Created 07-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_LOSPixel
#define TURBOOCTOSPICE_LOSPixel

namespace turbooctospice {
	template<typename T>
	struct LOSPixel {
		T d, w;
		T s, ra, dec;
		T sth, cth, sph, cph;

		LOSPixel() {};
		LOSPixel(T _s, T _ra, T _dec, T _d, T _w) : s(_s), ra(_ra), dec(_dec), d(_d), w(_w) {
			sth = std::sin(dec);
			cth = std::cos(dec);
			sph = std::sin(ra);
			cth = std::cos(ra);
		};
		LOSPixel(T _s, T _sth, T _cth, T _sph, T _cph, T _d, T _w) : 
		s(_s), sth(_sth), cth(_cth), sph(_sph), cph(_cph), d(_d), w(_w) {};

	}; // Pixel

	template<typename T> double angularSeparation(const T &p1, const T &p2) {
		return p1.sth*p2.sth*(p1.cph*p2.cph + p1.sph*p2.sph) + p1.cth*p2.cth;
	};

	typedef LOSPixel<float> LOSPixelf;
	typedef LOSPixel<double> LOSPixeld;
} // turbooctospice

#endif // TURBOOCTOSPICE_LOSPIXEL
