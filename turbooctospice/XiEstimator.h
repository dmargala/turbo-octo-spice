#ifndef TURBOOCTOSPICE_XI_ESTIMATOR
#define TURBOOCTOSPICE_XI_ESTIMATOR

#include <map>

#include <boost/progress.hpp>
#include <boost/timer.hpp>
#include <boost/thread/mutex.hpp>

#include "likely/likely.h"

#include "types.h"
#include "constants.h"

namespace turbooctospice {

    class XiEstimator {
    public:
        enum BinningCoordinateType {
            PolarCoordinates, CartesianCoordinates, ObservingCoordinates
        };
        XiEstimator(double scale, AbsTwoPointGridPtr grid,
            BinningCoordinateType type, std::vector<Forest> sightlines, SkyBinsIPtr skybins);
        void run(int nthreads);
        void save_results(std::string outfile);
        void save_subsamples(std::string outfile_base);
        void print_stats();
    private:
        void increment_progress();
        bool skybin_xi_task(int skybin_index);
        bool xi_finalize_task(int xi_bin_index);
        void accumulate_stats(unsigned long const &num_sightline_pairs,
            unsigned long const &num_sightline_pairs_used, unsigned long const &num_pixel_pairs,
            unsigned long const &num_pixel_pairs_used, unsigned long const &num_pixels);
        // private member variables
        // counters
        unsigned long num_sightline_pairs_, num_sightline_pairs_used_,
            num_pixels_, num_pixel_pairs_, num_pixel_pairs_used_;
        // inputs
        double max_ang_, cos_max_ang_;
        BinningCoordinateType coordinate_type_;
        AbsTwoPointGridPtr grid_;
        SkyBinsIPtr skybins_;
        std::vector<Forest> sightlines_;
        // results
        unsigned num_xi_bins_;
        std::map<int, std::vector<XiBin> > skybin_xis_;
        std::vector<XiBin> xi_;
        likely::CovarianceMatrixCPtr cov_matrix_;
        // helpers
        std::unique_ptr<boost::progress_display> show_progress_;
        boost::mutex show_progress_mutex_; // protects show_progress_
        boost::mutex pair_stats_mutex_; // protects counters
    };
} // turbooctospice

#endif // TURBOOCTOSPICE_XI_ESTIMATOR
