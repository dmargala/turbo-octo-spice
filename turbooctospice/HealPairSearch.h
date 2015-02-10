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
        typedef T PixelIterable;
        typedef typename PixelIterable::iterator PixelIterator;
        typedef std::map<int, std::vector<int> > BucketToPixels;
        HealPairSearch(PixelIterable quasars, HealMap map, double maxAng, bool verbose = false) : 
        _quasars(quasars), _map(map), _maxAng(maxAng), _verbose(verbose) {
            _cosmax = std::cos(_maxAng);
            initialize();
        };
        ~HealPairSearch() {};
        bool valid() const {
            if (_quasarIndexA < _quasars.size()) return true;
            return false;
        }
        void next() {
            ++_bPixelIt;
            if (_bPixelIt == _quasars[*_bQuasarIndexIt].pixels.end()) {
                ++_aPixelIt;
                if (_aPixelIt == _quasars[_quasarIndexA].pixels.end()) {
                    nextQuasarInBucket(); 
                    if(!valid()) return;
                    _newPair = true;
                    _aPixelIt = _quasars[_quasarIndexA].pixels.begin();
                }
                _bPixelIt = _quasars[*_bQuasarIndexIt].pixels.begin();
            }
        }
        template <class PairType> PairType get() {
            if (_newPair) {
                _cosij = PairType(_quasars[_quasarIndexA].pixels[0], _quasars[*_bQuasarIndexIt].pixels[0]).cosAngularSeparation();
                while(_cosij < _cosmax) {
                    next();
                    _cosij = PairType(_quasars[_quasarIndexA].pixels[0], _quasars[*_bQuasarIndexIt].pixels[0]).cosAngularSeparation();
                }
                _newPair = false;
            }
            return PairType(*_aPixelIt, *_bPixelIt);
        };

    private:
        void initialize() {
            for(int i = 0; i < _quasars.size(); ++i) {
                int index = _map.ang2pix(_quasars[i].p);
                if(_buckets.count(index) > 0 ) {
                    _buckets[index].push_back(i);
                }
                else {
                    _buckets[index] = std::vector<int>(1,i);
                }
            }
            _quasarIndexA = 0;
            resetNeighbors();
            nextBucketWithNeighbors();
            ////
            _bQuasarIndexIt = _buckets[*_neighboringBucketsIt].begin();
            nextQuasarInBucket();
            ////
            _aPixelIt = _quasars[_quasarIndexA].pixels.begin();
            _bPixelIt = _quasars[*_bQuasarIndexIt].pixels.begin();
            _newPair = true;
        }
        void nextQuasarInBucket(){
            // std::cerr << "\t\tnextQuasarInBucket (enter): ";
            // printState();
            ++_bQuasarIndexIt;
            if (_bQuasarIndexIt == _buckets[*_neighboringBucketsIt].end()) {
                nextBucketWithNeighbors();
                if(!valid()) return;
                _bQuasarIndexIt = _buckets[*_neighboringBucketsIt].begin();
            }
            if (*_bQuasarIndexIt <= _quasarIndexA) {
                nextQuasarInBucket();
            }
            // std::cerr << "\t\tnextQuasarInBucket (exit): ";
            // printState();
        }
        void nextBucketWithNeighbors() {
            // std::cerr << "\tnextBucketWithNeighbors (enter): ";
            // printState();
            ++_neighboringBucketsIt;
            if (_neighboringBucketsIt != _neighbors.end()) {
                if( !(_buckets.count(*_neighboringBucketsIt) > 0) ){
                    nextBucketWithNeighbors();
                }
            }
            else {
                nextQuasar();
                if (!valid()) return;
            }
            // std::cerr << "\tnextBucketWithNeighbors (exit): ";
            // printState();
        }
        void nextQuasar() {
            // std::cerr << "nextQuasar (enter): ";
            // printState();
            ++_quasarIndexA;
            if (!valid()) return;
            resetNeighbors();
            nextBucketWithNeighbors();
            // std::cerr << "nextQuasar (exit): ";
            // printState();
        }
        void resetNeighbors() {
            // Find neighboring buckets
            _map.query_disc_inclusive(_quasars[_quasarIndexA].p, _maxAng, _neighbors_rangeset);
            _neighbors_rangeset.toVector(_neighbors);
            _neighboringBucketsIt = _neighbors.begin();
        }
        void printState() {
            std::cerr << _quasarIndexA << " " 
                      << *_bQuasarIndexIt <<  " " 
                      << *_neighboringBucketsIt << " "
                      << std::endl;
        }
        /////////////////
        bool _verbose, _newPair;
        int _quasarIndexA;
        double _maxAng, _cosmax, _cosij;
        HealMap _map;
        PixelIterable _quasars;
        BucketToPixels _buckets;
        rangeset<int> _neighbors_rangeset;
        std::vector<int> _neighbors;
        std::vector<int>::iterator _neighboringBucketsIt;
        std::vector<int>::iterator _bQuasarIndexIt;
        std::vector<LOSPixelf>::iterator _aPixelIt, _bPixelIt;
    }; // HealPairSearch

    typedef HealPairSearch<std::vector<Quasarf> > HealSearch;
} // turbooctospice

#endif // TURBOOCTOSPICE_HEAL_PAIR_SEARCH
