// Created 28-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>

#include "MockSpectrum.h"

#include "boost/foreach.hpp"
#include "boost/smart_ptr.hpp"
#include "boost/format.hpp"
#include "boost/algorithm/string.hpp"

#include <CCfits/CCfits>

#include <stdexcept>
#include <algorithm>
#include <string>

namespace local = turbooctospice;

const double PI = std::atan(1.0)*4;
const double DEG2RAD = PI/180.0;

local::MockSpectrum::MockSpectrum(std::string target, int id, bool verbose) :
_target(target), _id(id) {
    loadTarget(verbose);
}

void local::MockSpectrum::loadTarget(bool verbose) {
    try {
        // Lookup file name for this target
        std::string rawFilename = getMockFilename(_target);
        if(verbose) {
            std::cout << "Reading raw mock file " << rawFilename << std::endl;
        }
        // Read fits file
        std::unique_ptr<CCfits::FITS> rawMockFile(new CCfits::FITS(rawFilename, CCfits::Read));
        // Read header keywords
        CCfits::PHDU& header = rawMockFile->pHDU();
        header.readKey("m_z", _z);
        header.readKey("m_ra", _ra);
        header.readKey("m_dec", _dec);
        header.readKey("coeff0", _coeff0);
        header.readKey("coeff1", _coeff1);
        // Read mock extension
        CCfits::ExtHDU& table = rawMockFile->extension(1);
        table.column("f").read(_frac, 1, table.rows());
    }
    catch (CCfits::FitsException& e) {
        std::cerr << "CCfits Exception Thrown :" << e.message();
    }
}

local::Forest local::MockSpectrum::getForest(
    int ncombine, float forestlo, float foresthi, float speclo) {
    // init forest pixels
    Forest forest(_ra*DEG2RAD, _dec*DEG2RAD, _id);
    // double theta((90.0-_dec)*DEG2RAD), phi(_ra*DEG2RAD);
    // forest.phi = phi;
    // forest.theta = theta;
    // forest.sdec = std::sin(_dec*DEG2RAD);
    // forest.cdec = std::cos(_dec*DEG2RAD);
    // forest.sph = std::sin(phi);
    // forest.cph = std::cos(phi);
    // Max and min wavelength defined by lyman alpha forest range and spec cutoff
    float minLogLambda(std::log10(std::max(forestlo*(1+_z), speclo)));
    float maxLogLambda(std::log10(foresthi*(1+_z)));
    // Iterate over raw mock pixels
    float logLambda, frac, weight, distance;
    for(int i = 0; i < _frac.size()-ncombine+1; i+=ncombine) {
        logLambda = 0; frac = 0; weight = 0; distance = 0;
        for(int j = 0; j < ncombine; ++j) {
            logLambda += _coeff0+_coeff1*(i+j);
            frac += _frac[i+j];
            weight += 1.0;
            distance += 0;
        }
        logLambda /= ncombine;
        // If below forest, skip to next pixel
        if (logLambda < minLogLambda) continue;
        // If beyond forest, we're done
        if (logLambda > maxLogLambda) break;
        // Save the pixel
        forest.pixels.push_back( {frac/ncombine, logLambda, weight/ncombine, distance/ncombine} );
    }
    return forest;
}


std::string local::getMockFilename(std::string target) {
    const char *fromEnv(std::getenv("BOSS_ROOT"));
    if(0 == fromEnv) {
        throw std::runtime_error("Environment variable BOSS_ROOT is not set.");
    }
    std::string BOSS_ROOT = fromEnv;
    fromEnv = std::getenv("MOCK_VERSION");
    if(0 == fromEnv) {
        throw std::runtime_error("Environment variable MOCK_VERSION is not set.");
    }
    std::string MOCK_VERSION = fromEnv;

    boost::format filenameFormat("%s/%s/rawlite/%s/mockrawShort-%s.fits");
    std::vector<std::string> strs;
    boost::split(strs, target, boost::is_any_of("-"));
    return boost::str(filenameFormat % BOSS_ROOT % MOCK_VERSION % strs[0] % target);
}
