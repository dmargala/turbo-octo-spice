// Created 07-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_PAIR
#define TURBOOCTOSPICE_PAIR

#include "utility"

#include "boost/coroutine/coroutine.hpp"

namespace turbooctospice {
	template <typename T>
	class Pair : public std::pair<T,T> {
	public:
		typedef std::pair<T,T> Base;
		typedef boost::coroutines::coroutine<Pair> PairGenerator;
		Pair() : Base() {};
		Pair(const Pair &pair) : Base(pair) {};
		Pair(T first, T second) : Base(first,second) {};
		double separation(){
			double dx = Base::first.x - Base::second.x;
            double dy = Base::first.y - Base::second.y;
            double dz = Base::first.z - Base::second.z;
            return std::sqrt(dx*dx+dy*dy+dz*dz);
		}
		double separationSq(){
			double dx = Base::first.x - Base::second.x;
            double dy = Base::first.y - Base::second.y;
            double dz = Base::first.z - Base::second.z;
            return dx*dx+dy*dy+dz*dz;
		}
		double weight(){
            return Base::first.w*Base::second.w;
		}
		double product(){
			return Base::first.d*Base::second.d;
		}
	private:
	}; // Pair

	typedef Pair<Pixel> XYZPixelPair;

	template <typename T>
	class AngPair : public std::pair<T,T> {
	public:
		typedef std::pair<T,T> Base;
		typedef boost::coroutines::coroutine<AngPair> PairGenerator;

		AngPair() : Base() {};
		AngPair(const AngPair &pair) : Base(pair) { 
			_cos12 = angularSeparation<T>(pair.first, pair.second);
		};
		AngPair(T first, T second) : Base(first,second) {
			_cos12 = angularSeparation<T>(first, second);
		};
		AngPair(T first, T second, double cos12) : Base(first,second), _cos12(cos12) {};
		double separation(){
			double s1 = Base::first.s;
			double s2 = Base::second.s;
			return std::sqrt(s1*s1 + s2*s2 - 2*s1*s2*_cos12);
		}
		double separationSq(){
			double s1 = Base::first.s;
			double s2 = Base::second.s;
			return s1*s1 + s2*s2 - 2*s1*s2*_cos12;
		}
		double cosAngularSeparation(){
			return _cos12;
		}
		double weight(){
            return Base::first.w*Base::second.w;
		}
		double product(){
			return Base::first.d*Base::second.d;
		}
	private:
		double _cos12;
	}; // AngPair

	typedef AngPair<LOSPixelf> AngPixelPair;

} // turbooctospice

#endif // TURBOOCTOSPICE_PAIR
