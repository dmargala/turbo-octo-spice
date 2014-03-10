// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_BIN_XYZ_PAIR
#define TURBOOCTOSPICE_BIN_XYZ_PAIR

#include "types.h"

#include "likely/likely.h"

#include <iostream>

namespace lk = likely;

namespace turbooctospice {

	template <class T>
	class BinXYZPair {
	public:
        typedef T PairType;
        BinXYZPair(double minValue, double maxValue, int nBins, bool countPair = false): 
        _minValue(minValue), _maxValue(maxValue), _nBins(nBins), _countPair(countPair) {
            _binWidth = (maxValue - minValue)/nBins;
            _minValueSq = _minValue*_minValue;
            _maxValueSq = _maxValue*_maxValue;
        };
        ~BinXYZPair() {};
        //template <class PairType> void binPair(typename PairType &pair, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
		void binPair(PairType &pair, std::vector<double> &dsum, std::vector<double> &wsum, long &nused) const {
            // Check that separation is within range of interest
            double rsq = pair.separationSq();
            if(rsq < _minValueSq || rsq >= _maxValueSq) return;
            // Accumulate pair
            int index = getBinIndex(std::sqrt(rsq));
            double wgt = pair.weight();
            wsum[index] += wgt;
            if(_countPair) {
                dsum[index] = wsum[index]*wsum[index];
            }
            else {
                dsum[index] += wgt*pair.product();
            }
            nused++;
		};
        int getNBins() const { return _nBins; };
        int getBinIndex(double value) const {
            return std::floor((value - _minValue)/_binWidth);
        }
        double getBinCenter(int index) const {
            return _minValue + (index+0.5)*_binWidth;
        }
	private:
        double _minValue, _maxValue, _nBins, _binWidth, _minValueSq, _maxValueSq;
        bool _countPair;
	}; // BinXYZPair

    typedef BinXYZPair<XYZPixelPair> BinXYZ;
    typedef BinXYZPair<AngPixelPair> BinAng;

} // turbooctospice

#endif // TURBOOCTOSPICE_BIN_XYZ_PAIR
