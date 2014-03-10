// Created 07-Mar-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#ifndef TURBOOCTOSPICE_HEAL_PAIR_SEARCH
#define TURBOOCTOSPICE_HEAL_PAIR_SEARCH

#include "boost/shared_ptr.hpp"

#include "boost/progress.hpp"

#include "healpix_map.h" // Healpix_Map, rangeset, pointing
#include <cmath>

namespace turbooctospice {

	typedef Healpix_Map<double> HealMap;
    typedef boost::shared_ptr<HealMap> HealMapPtr;

    template <class PixelType>
    struct Quasar {
        std::vector<PixelType> pixels;
        pointing p;
    };

    typedef Quasar<LOSPixelf> Quasarf;
    typedef std::vector<Quasarf> Quasars;


	template <class T>
	class HealPairSearch {
	public:
		typedef std::map<int, std::vector<int> > BucketToPixels;
		typedef T PixelIterable;
		HealPairSearch(PixelIterable quasars, HealMap map, double maxAng, bool verbose = false) : 
		_quasars(quasars), _map(map), _maxAng(maxAng), _verbose(verbose) {
			_cosmax = std::cos(_maxAng);
		};
		~HealPairSearch() {};
        template <class PairGenerator, class PairType> void findPairs(typename PairGenerator::caller_type& yield) const {
        	if (_verbose) std::cout << "Entering auto-correlation generator ..." << std::endl;
		    rangeset<int> neighbors_rangeset;
		    std::vector<int> neighbors;

			BucketToPixels _buckets;
		    sortQuasars(_buckets);

    		if (_verbose) std::cout << "We have " << _buckets.size() << " buckets w/ data" << std::endl;

    		boost::progress_display t(_quasars.size());
	        // Loop over all quasars in this bucket
	        for(int quasarIndexI = 0; quasarIndexI < _quasars.size(); ++quasarIndexI) {
	            // Find neighboring buckets
	            _map.query_disc_inclusive(_quasars[quasarIndexI].p, _maxAng, neighbors_rangeset);
	            neighbors_rangeset.toVector(neighbors);
	            // Loop over neighboring buckets
	            for(int neighborBucketIndex : neighbors) {
	                // Loop over all quasars in neighboring bucket
	                for(int quasarIndexJ : _buckets[neighborBucketIndex]) {
	                    // Only use unique pairs of quasars
	                    if(quasarIndexJ <= quasarIndexI) continue;
	                    // kinda kuldgey...
	                    double cosij = PairType(_quasars[quasarIndexI].pixels[0], _quasars[quasarIndexJ].pixels[0]).cosAngularSeparation();
	                    if(cosij < _cosmax) continue;
	                    for(int i = 0; i < _quasars[quasarIndexI].pixels.size(); ++i) {
	                        for(int j = 0; j < _quasars[quasarIndexJ].pixels.size(); ++j) {
                            	yield(PairType(_quasars[quasarIndexI].pixels[i],_quasars[quasarIndexJ].pixels[j],cosij));
	                        }
	                    }
	                }
	            }
	            t+=1;
	        }
		    if (_verbose) std::cout << "Exiting auto-correlation generator ..." << std::endl;
        }
	private:
		void sortQuasars(BucketToPixels &_buckets) const {
			for(int i = 0; i < _quasars.size(); ++i) {
				int index = _map.ang2pix(_quasars[i].p);
				if(_buckets.count(index) > 0 ) {
					_buckets[index].push_back(i);
				}
				else {
					_buckets[index] = std::vector<int>(1,i);
				}
			}
		}
		bool _verbose;
		double _maxAng, _cosmax;
		HealMap _map;
		PixelIterable _quasars;
	}; // HealPairSearch

	typedef HealPairSearch<std::vector<Quasarf>> HealSearch;
} // turbooctospice

#endif // TURBOOCTOSPICE_HEAL_PAIR_SEARCH
