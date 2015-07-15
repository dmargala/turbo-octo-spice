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
        bool healxi_task(int id);
        bool xi_finalize_task(int id);
        bool cov_task(int a, int b);
        void accumulate_stats(unsigned long const &num_sightline_pairs,
            unsigned long const &num_sightline_pairs_used, unsigned long const &num_pixel_pairs,
            unsigned long const &num_pixel_pairs_used, unsigned long const &num_pixels);
        unsigned num_xi_bins_;
        unsigned long num_pixels_, num_sightlines_,
            num_sightline_pairs_, num_sightline_pairs_used_,
            num_pixel_pairs_, num_pixel_pairs_used_;
        double max_ang_, cos_max_ang_;
        BinningCoordinateType coordinate_type_;
        SkyBinsIPtr skybins_;
        AbsTwoPointGridPtr grid_;
        std::vector<Forest> sightlines_;
        std::map<int, std::vector<XiBin> > healxis_;
        std::vector<XiBin> xi_;
        std::vector<std::vector<double> > cov_;
        std::unique_ptr<boost::progress_display> show_progress_;
        boost::mutex show_progress_mutex_; // protects show_progress_
        boost::mutex pair_stats_mutex_;
        boost::mutex debug_mutex_;
        likely::CovarianceMatrixCPtr cov_matrix_;
    };
} // turbooctospice

#endif // TURBOOCTOSPICE_XI_ESTIMATOR
