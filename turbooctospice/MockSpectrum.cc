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

local::MockSpectrum::MockSpectrum(std::string target, bool verbose) : 
_target(target) {
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
        std::auto_ptr<CCfits::FITS> rawMockFile(new CCfits::FITS(rawFilename,CCfits::Read));
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

std::vector<local::QuasarPixel> local::MockSpectrum::getTrimmedSpectrum(
int ncombine, float forestlo, float foresthi, float speclo) {
    // init quasar pixels
    std::vector<QuasarPixel> pixels;
    // Max and min wavelength defined by lyman alpha forest range and spec cutoff
    float minlam(std::max(forestlo*(1+_z), speclo));
    float maxlam(foresthi*(1+_z));
    // Iterate over raw mock pixels
    float lam;
    QuasarPixel qp;
    for(int i = 0; i < _frac.size(); ++i) {
        // Calculate current pixel's wavelength
        lam = std::pow(10,_coeff0+i*_coeff1);
        // If below forest, skip to next pixel
        if (lam < minlam) continue;
        // If beyond forest, we're done
        if (lam > maxlam) break;
        // Save the pixel
        qp.frac = _frac[i];
        qp.lam = lam;
        qp.wgt = 1.0;
        pixels.push_back(qp);
    }
    return pixels;
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