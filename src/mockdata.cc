// Created 24-Feb-2014 by Daniel Margala (University of California, Irvine) <dmargala@uci.edu>
// Read BOSS mock data 

#include "boost/program_options.hpp"
#include "boost/foreach.hpp"
#include "boost/format.hpp"
#include "boost/smart_ptr.hpp"

#include "bosslya/bosslya.h"
#include "cosmo/cosmo.h"
#include "likely/likely.h"

#include <iostream>
#include <fstream>
#include <algorithm>

#include <CCfits/CCfits>
#include "healpix_map.h"


namespace po = boost::program_options;
namespace lk = likely;
namespace lya = bosslya;

struct QuasarPixel {
    float flux, lam, wgt, dist;
};

struct QuasarSpectrum {
    float z, ra, dec, sindec, cosdec, sinra, cosra;
    pointing p;
    std::vector<QuasarPixel> pixels;
};


int main(int argc, char **argv) {

    std::string infile;
    po::options_description cli("Read mock data");
    cli.add_options()
        ("help,h", "Prints this info and exits.")
        ("verbose", "Prints additional information.")
        ("input,i", po::value<std::string>(&infile)->default_value(""),
            "Filename to read from")
        ;
    // do the command line parsing now
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, cli), vm);
        po::notify(vm);
    }
    catch(std::exception const &e) {
        std::cerr << "Unable to parse command line options: " << e.what() << std::endl;
        return -1;
    }
    if(vm.count("help")) {
        std::cout << cli << std::endl;
        return 1;
    }
    bool verbose(vm.count("verbose"));

    long nTargets;
    lya::TargetSet targetSet;

    // Read targetlist file with redshifts
    nTargets = lya::readTargetSet(targetSet,infile);
    if(verbose) {
        std::cout << "Will read " << nTargets << " targets." << std::endl;
    }

    std::string rawFilename;

    boost::scoped_ptr<CCfits::FITS> pinfile;

    std::string baseDir = "/Users/daniel/data/boss";
    std::string mockDir = "M3_0_0/rawlite";
    boost::format filenameFormat("%s/%s/%04d/mockrawShort-%04d-%5d-%04d.fits");

    std::vector<QuasarSpectrum> quasars;

    // Loop over the available targets in the input file.        
    BOOST_FOREACH(lya::Target target, targetSet) {
        try {
            rawFilename = (filenameFormat % baseDir % mockDir % target.getPlate() % target.getPlate() % target.getMJD() % target.getFiber()).str();

            pinfile.reset(new CCfits::FITS(rawFilename,CCfits::Read));

            CCfits::PHDU& header = pinfile->pHDU();

            if(verbose) {
                std::cout << "Reading raw mock file " << rawFilename << std::endl;
            }

            QuasarSpectrum quasar;

            header.readKey("m_z", quasar.z);
            header.readKey("m_ra", quasar.ra);
            header.readKey("m_dec", quasar.dec);
            double coeff0, coeff1;
            header.readKey("coeff0", coeff0);
            header.readKey("coeff1", coeff1);

            boost::format outFormat("%.4f %.4f %.4f %.4f %.4f");
            std::cout << outFormat % quasar.z % quasar.ra % quasar.dec % coeff0 % coeff1 << std::endl;

            CCfits::ExtHDU& table = pinfile->extension(1);
            std::vector<float> f;

            table.column("f").read(f, 1, table.rows());

            float lam;
            float forestlo(1040), foresthi(1200), lya(1216), speclo(3650);

            float minlam = std::max(forestlo*(1+quasar.z), speclo);
            float maxlam = foresthi*(1+quasar.z);

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

            std::cout << f.size() << " " << quasar.pixels.size() << std::endl;
        }
        catch (CCfits::FitsException&) 
        {
          std::cerr << " Fits Exception Thrown by test function \n";       
        }
    }

    // Read the input file
    if(0 == infile.length()) {
        std::cerr << "Missing infile parameter." << std::endl;
        return -2;
    }

    return 0;
}