#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <map>

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/thread/mutex.hpp>

#include "cosmo/cosmo.h"

#include "types.h"
#include "constants.h"

namespace turbooctospice {

    class XiEstimator {
    public:
        enum BinningCoordinateType { PolarCoordinates, CartesianCoordinates, ObservingCoordinates };
        XiEstimator(int order, std::string infile, cosmo::AbsHomogeneousUniversePtr cosmology,
            AbsTwoPointGridPtr grid, BinningCoordinateType type, bool skip_ngc=false, bool skip_sgc=false);
        void run(int nthreads);
        void save_results(std::string outfile);
    private:
        // void task_finalize();
        void increment_progress();
        bool healxi_task(int id);
        bool xi_finalize_task(int id);
        bool cov_task(int a, int b);

        int num_xi_bins, num_mu_bins, num_rperp_bins, num_sep_bins, num_z_bins;
        float r_min, r_max, r_spacing, rsq_min, rsq_max,
            rperp_min, rperp_max, rperp_spacing, sep_min, sep_spacing,
            z_min, z_max, z_spacing;
        double max_ang_, cos_max_ang_;
        BinningCoordinateType coordinate_type_;
        HealpixBinsI healbins_;
        AbsTwoPointGridPtr grid_;
        std::vector<Forest> sightlines_;
        std::map<int, std::vector<XiBin> > healxis_;
        std::vector<XiBin> xi_;
        std::vector<std::vector<double> > cov_;
        boost::progress_display show_progress_;
        boost::mutex show_progress_mutex_; // protects show_progress_

    };
} // turbooctospice

#endif // TURBOOCTOSPICE_XI_ESTIMATOR
