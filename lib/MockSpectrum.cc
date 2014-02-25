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

local::MockSpectrum::MockSpectrum(float z, float ra, float dec) : 
_z(z), _ra(ra), _dec(dec) {

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

void local::readMockTargets(std::vector<std::string> &targetlist, std::vector<MockSpectrum> &quasars, bool const verbose=false) {
    boost::scoped_ptr<CCfits::FITS> rawMockFile;
    // Loop over the available targets in the input file.        
    BOOST_FOREACH(std::string target, targetlist) {
        try {
            std::string rawFilename = getMockFilename(target);
            if(verbose) {
                std::cout << "Reading raw mock file " << rawFilename << std::endl;
            }
            rawMockFile.reset(new CCfits::FITS(rawFilename,CCfits::Read));

            CCfits::PHDU& header = rawMockFile->pHDU();

            float z, ra, dec;
            header.readKey("m_z", z);
            header.readKey("m_ra", ra);
            header.readKey("m_dec", dec);
            MockSpectrum quasar(z, ra, dec);

            float coeff0, coeff1;
            header.readKey("coeff0", coeff0);
            header.readKey("coeff1", coeff1);

            CCfits::ExtHDU& table = rawMockFile->extension(1);
            std::vector<float> f;

            table.column("f").read(f, 1, table.rows());

            float lam;
            float forestlo(1040), foresthi(1200), lya(1216), speclo(3650);

            float minlam(std::max(forestlo*(1+z), speclo));
            float maxlam(foresthi*(1+z));

            QuasarPixel qp;
            for(int i = 0; i < f.size(); ++i) {
                lam = std::pow(10,coeff0+i*coeff1);
                if (lam < minlam) continue;
                if (lam > maxlam) break;
                qp.flux = f[i];
                qp.lam = lam;
                qp.wgt = 1.0;
                quasar.pixels.push_back(qp);
            }
            quasars.push_back(quasar);
        }
        catch (CCfits::FitsException& e) {
            std::cerr << "CCfits Exception Thrown :" << e.message();       
        }
    }
}
